//#####################################################################
// Copyright 2003, 2003, Ron Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFT_2D
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms/FFT_2D.h>
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
template<class T> FFT_2D<T>::
FFT_2D(const GRID<TV>& grid_input)
    :grid(grid_input),plan_u_to_u_hat(0),plan_u_hat_to_u(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FFT_2D<T>::
~FFT_2D()
{
    FFTW<T,2>::Destroy_Plan(plan_u_to_u_hat);
    FFTW<T,2>::Destroy_Plan(plan_u_hat_to_u);
}
//#####################################################################
// Function Transform
//#####################################################################
template<class T> void FFT_2D<T>::
Transform(const ARRAY<T,VECTOR<int,2> >& u,ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat) const
{
    if(plan_u_to_u_hat_counts!=grid.Counts()){
        plan_u_to_u_hat_counts=grid.Counts();
        FFTW<T,2>::Destroy_Plan(plan_u_to_u_hat);
        plan_u_to_u_hat=FFTW<T,2>::Plan_R2C(grid.Counts(),u,u_hat);}
    FFTW<T,2>::Execute_R2C(plan_u_to_u_hat,u,u_hat);
}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class T> void FFT_2D<T>::
Inverse_Transform(ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat,ARRAY<T,VECTOR<int,2> >& u,bool normalize,bool preserve_u_hat) const
{
    if(plan_u_hat_to_u_counts!=grid.Counts()){
        plan_u_hat_to_u_counts=grid.Counts();
        FFTW<T,2>::Destroy_Plan(plan_u_hat_to_u);
        plan_u_hat_to_u=FFTW<T,2>::Plan_C2R(grid.Counts(),u_hat,u);}
    if(preserve_u_hat){
        u_hat_copy=u_hat;
        FFTW<T,2>::Execute_C2R(plan_u_hat_to_u,u_hat_copy,u);}
    else FFTW<T,2>::Execute_C2R(plan_u_hat_to_u,u_hat,u);
    if(normalize) u*=(T)1/grid.Counts().Product();
}
//#####################################################################

#else // NR-specific functions **********************************************************************************************************************************************************************************
#include <PhysBAM_Tools/Fourier_Transforms/FFT.h>
//#####################################################################
// Constructor
//#####################################################################
template<class T> FFT_2D<T>::
FFT_2D(const GRID<TV>& grid_input)
    :grid(grid_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FFT_2D<T>::
~FFT_2D()
{}
//#####################################################################
// Function Transform
//#####################################################################
template<class T> void FFT_2D<T>::
Transform(const ARRAY<T,VECTOR<int,2> >& u,ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat) const
{
    ARRAY<int> dim(grid.counts);
    data.Resize(2*grid.counts.x*grid.counts.y,false,false);
    for(int i=1,k=0;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){data(++k)=(float)u(i,j);data(++k)=0;}
    NR_fourn(-1,dim,data);
    for(int i=0,k=0;i<=grid.counts.x-1;i++) {for(int j=0;j<=grid.counts.y/2;j++){u_hat(i,j).re=data(++k);u_hat(i,j).im=data(++k);} k+=grid.counts.y-2;}
}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class T> void FFT_2D<T>::
Inverse_Transform(ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat,ARRAY<T,VECTOR<int,2> >& u,bool normalize,bool preserve_u_hat) const
{
    ARRAY<int> dim(grid.counts);
    data.Resize(2*grid.counts.x*grid.counts.y,false,false);
    int k=0;
    for(int i=0;i<=grid.counts.x-1;i++){int negi=(i==0)?0:grid.counts.x-i;
        for(int j=0;j<=grid.counts.y/2;j++){data(++k)=(float)u_hat(i,j).re;data(++k)=(float)u_hat(i,j).im;}
        for(int j=grid.counts.y/2+1;j<=grid.counts.y-1;j++){data(++k)=(float)u_hat(negi,grid.counts.y-j).re;data(++k)=-(float)u_hat(negi,grid.counts.y-j).im;}}
    NR_fourn(+1,dim,data);
    if(normalize){T coefficient=T(1)/(grid.counts.x*grid.counts.y);for(int k=0,i=1;i<=grid.counts.x;i++)for(int j=1;j<=grid.counts.y;j++){u(i,j)=coefficient*(T)data(++k);++k;}}
    else for(int k=0,i=1;i<=grid.counts.x;i++)for(int j=1;j<=grid.counts.y;j++){u(i,j)=(T)data(++k);++k;}
}
//#####################################################################
#endif // ********************************************************************************************************************************************************************************************************

//#####################################################################
// Function Enforce_Real_Valued_Symmetry
//#####################################################################
// enforce symmetry so that the inverse transform is a real valued function
template<class T> void FFT_2D<T>::
Enforce_Real_Valued_Symmetry(ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat) const
{
    // imaginary part of the constant and cosine only terms are identically zero
    u_hat(0,0).im=u_hat(grid.counts.x/2,0).im=u_hat(0,grid.counts.y/2).im=u_hat(grid.counts.x/2,grid.counts.y/2).im=0;
    // enforce appropriate complex conjugates
    for(int i=1;i<=grid.counts.x/2-1;i++){
        u_hat(grid.counts.x-i,0)=u_hat(i,0).Conjugated(); // reflect y=0 line
        u_hat(grid.counts.x-i,grid.counts.y/2)=u_hat(i,grid.counts.y/2).Conjugated();} // reflect y=grid.counts.y/2 line
}
//#####################################################################
// Function First_Derivatives
//#####################################################################
template<class T> void FFT_2D<T>::
First_Derivatives(const ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat,ARRAY<COMPLEX<T> ,VECTOR<int,2> >& ux_hat,ARRAY<COMPLEX<T> ,VECTOR<int,2> >& uy_hat) const
{
    VECTOR<T,2> coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    for(int i=0;i<=grid.counts.x-1;i++){T k1=coefficients.x*(i<=grid.counts.x/2?i:i-grid.counts.x);
        for(int j=0;j<=grid.counts.y/2;j++){T k2=coefficients.y*j;
            ux_hat(i,j)=uy_hat(i,j)=u_hat(i,j).Rotated_Counter_Clockwise_90();
            ux_hat(i,j)*=k1;uy_hat(i,j)*=k2;}}
    Enforce_Real_Valued_Symmetry(ux_hat);Enforce_Real_Valued_Symmetry(uy_hat);
}
//#####################################################################
// Make_Divergence_Free
//#####################################################################
// the constant fields are unchanged
template<class T> void FFT_2D<T>::
Make_Divergence_Free(ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat,ARRAY<COMPLEX<T> ,VECTOR<int,2> >& v_hat) const
{
    VECTOR<T,2> coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    // u=0 on the bottom
    for(int i=1;i<=grid.counts.x/2;i++) u_hat(i,0)=COMPLEX<T>(0,0);
    // area
    for(int i=0;i<=grid.counts.x-1;i++){T k1=coefficients.x*(i<=grid.counts.x/2?i:i-grid.counts.x);
        for(int j=1;j<=grid.counts.y/2;j++){T k2=coefficients.y*j,one_over_k_squared=1/(sqr(k1)+sqr(k2));
            COMPLEX<T> correction=(k1*u_hat(i,j)+k2*v_hat(i,j))*one_over_k_squared;
            u_hat(i,j)-=correction*k1;v_hat(i,j)-=correction*k2;}}
    Enforce_Real_Valued_Symmetry(u_hat);Enforce_Real_Valued_Symmetry(v_hat);
}
//#####################################################################
// Function Filter_High_Frequencies
//#####################################################################
// Lanczos filter - doesn't change (0,0) frequency
template<class T> void FFT_2D<T>::
Filter_High_Frequencies(ARRAY<COMPLEX<T> ,VECTOR<int,2> >& u_hat,T scale) const
{
    T coefficient=2*T(pi)/sqrt((T)(sqr(grid.counts.x)+sqr(grid.counts.y)));
    for(int i=1;i<=grid.counts.x/2;i++){ // skip the (0,0) case
        T temp=scale*coefficient*i,damping=sinc(temp);
        u_hat(i,0)*=damping;}
    for(int i=0;i<=grid.counts.x-1;i++){T i_frequency=T(i<=grid.counts.x/2 ? i:i-grid.counts.x);
        for(int j=1;j<=grid.counts.y/2;j++){T j_frequency=T(j);
            T temp=scale*coefficient*sqrt(sqr(i_frequency)+sqr(j_frequency)),damping=sinc(temp);
            u_hat(i,j)*=damping;}}
    Enforce_Real_Valued_Symmetry(u_hat);
}
//#####################################################################
template class FFT_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FFT_2D<double>;
#endif
