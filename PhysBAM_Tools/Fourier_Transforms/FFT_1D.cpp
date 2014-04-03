//#####################################################################
// Copyright 2002-2006, Ron Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFT_1D
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms/FFT_1D.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFTW.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Functions.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
using namespace PhysBAM;

#ifdef USE_FFTW // FFTW-specific functions *******************************************************************************************************************************************************************
//#####################################################################
// Constructor
//#####################################################################
template<class T> FFT_1D<T>::
FFT_1D(const GRID<TV>& grid_input)
    :grid(grid_input),plan_u_to_u_hat(0),plan_u_hat_to_u(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FFT_1D<T>::
~FFT_1D()
{
    FFTW<T,1>::Destroy_Plan(plan_u_to_u_hat);
    FFTW<T,1>::Destroy_Plan(plan_u_hat_to_u);
}
//#####################################################################
// Function Transform
//#####################################################################
template<class T> void FFT_1D<T>::
Transform(const ARRAY<T,VECTOR<int,1> >& u,ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat) const
{
    if(plan_u_to_u_hat_counts!=grid.Counts()){
        plan_u_to_u_hat_counts=grid.Counts();
        FFTW<T,1>::Destroy_Plan(plan_u_to_u_hat);
        plan_u_to_u_hat=FFTW<T,1>::Plan_R2C(grid.Counts(),u,u_hat);}
    FFTW<T,1>::Execute_R2C(plan_u_to_u_hat,u,u_hat);
}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class T> void FFT_1D<T>::
Inverse_Transform(ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat,ARRAY<T,VECTOR<int,1> >& u,bool normalize,bool preserve_u_hat) const
{
    if(plan_u_hat_to_u_counts!=grid.Counts()){
        plan_u_hat_to_u_counts=grid.Counts();
        FFTW<T,1>::Destroy_Plan(plan_u_hat_to_u);
        plan_u_hat_to_u=FFTW<T,1>::Plan_C2R(grid.Counts(),u_hat,u);}
    if(preserve_u_hat){
        u_hat_copy=u_hat;
        FFTW<T,1>::Execute_C2R(plan_u_hat_to_u,u_hat_copy,u);}
    else FFTW<T,1>::Execute_C2R(plan_u_hat_to_u,u_hat,u);
    if(normalize) u*=(T)1/grid.Counts().Product();
}
//#####################################################################

#else // NR-specific functions **********************************************************************************************************************************************************************************
#include <PhysBAM_Tools/Fourier_Transforms/FFT.h>

//#####################################################################
// Constructor
//#####################################################################
template<class T> FFT_1D<T>::
FFT_1D(const GRID<TV>& grid_input)
    :grid(grid_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FFT_1D<T>::
~FFT_1D()
{}
//#####################################################################
// Function Transform
//#####################################################################
template<class T> void FFT_1D<T>::
Transform(const ARRAY<T,VECTOR<int,1> >& u,ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat) const
{
    ARRAY<int> dim(grid.counts);
    data.Resize(2*grid.counts.x,false,false);
    for(int i=1,k=0;i<=grid.counts.x;i++){data(++k)=(float)u(i);data(++k)=0;}
    NR_fourn(-1,dim,data);
    for(int i=0,k=0;i<=grid.counts.x/2;i++){u_hat(i).re=data(++k);u_hat(i).im=data(++k);}
}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class T> void FFT_1D<T>::
Inverse_Transform(ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat,ARRAY<T,VECTOR<int,1> >& u,bool normalize,bool preserve_u_hat) const
{
    ARRAY<int> dim(grid.counts);
    data.Resize(2*grid.counts.x,false,false);
    int k=0;
    for(int i=0;i<=grid.counts.x/2;i++){data(++k)=(float)u_hat(i).re;data(++k)=(float)u_hat(i).im;}
    for(int i=grid.counts.x/2+1;i<=grid.counts.x-1;i++){data(++k)=(float)u_hat(grid.counts.x-i).re;data(++k)=-(float)u_hat(grid.counts.x-i).im;}
    NR_fourn(+1,dim,data);
    if(normalize){T coefficient=(T)1/grid.counts.x;k=0;for(int i=1,k=0;i<=grid.counts.x;i++){u(i)=coefficient*data(++k);++k;}}
    else for(int i=1,k=0;i<=grid.counts.x;i++){u(i)=data(++k);++k;}
}
//#####################################################################
#endif

//#####################################################################
// Function Enforce_Real_Valued_Symmetry
//#####################################################################
// enforce symmetry so that the inverse transform is a real valued function
template<class T> void FFT_1D<T>::
Enforce_Real_Valued_Symmetry(ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat) const
{
    // imaginary part of the constant and cosine only terms are identically zero
    u_hat(0).im=u_hat(grid.counts.x/2).im=0;
}
//#####################################################################
// Function First_Derivatives
//#####################################################################
template<class T> void FFT_1D<T>::
First_Derivatives(const ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat,ARRAY<COMPLEX<T> ,VECTOR<int,1> >& ux_hat) const
{
    VECTOR<T,1> coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    for(int i=0;i<=grid.counts.x/2;i++){T k=coefficients.x*i;
        ux_hat(i)=k*u_hat(i).Rotated_Counter_Clockwise_90();}
    Enforce_Real_Valued_Symmetry(ux_hat);
}
//#####################################################################
// Make_Divergence_Free
//#####################################################################
// only the constant field is nonzero in 1D
template<class T> void FFT_1D<T>::
Make_Divergence_Free(ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat) const
{
    for(int i=1;i<=grid.counts.x/2;i++) u_hat(i)=COMPLEX<T>(0,0);
    Enforce_Real_Valued_Symmetry(u_hat);
}
//#####################################################################
// Function Filter_High_Frequencies
//#####################################################################
// Lanczos filter - doesn't change (0,0) frequency
template<class T> void FFT_1D<T>::
Filter_High_Frequencies(ARRAY<COMPLEX<T> ,VECTOR<int,1> >& u_hat,T scale) const
{
    T coefficient=2*(T)pi/grid.counts.x;
    for(int i=1;i<=grid.counts.x/2;i++){
        T temp=scale*coefficient*i,damping=sinc(temp);
        u_hat(i)*=damping;}
    Enforce_Real_Valued_Symmetry(u_hat);
}
//#####################################################################
template class FFT_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FFT_1D<double>;
#endif
