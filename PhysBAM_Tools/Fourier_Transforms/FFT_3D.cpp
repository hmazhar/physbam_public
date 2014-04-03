//#####################################################################
// Copyright 2002-2006, Ron Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFT_3D
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms/FFT_3D.h>
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
template<class T> FFT_3D<T>::
FFT_3D(const GRID<TV>& grid_input)
    :grid(grid_input),plan_u_to_u_hat(0),plan_u_hat_to_u(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FFT_3D<T>::
~FFT_3D()
{
    FFTW<T,3>::Destroy_Plan(plan_u_to_u_hat);
    FFTW<T,3>::Destroy_Plan(plan_u_hat_to_u);
}
//#####################################################################
// Function Transform
//#####################################################################
template<class T> void FFT_3D<T>::
Transform(const ARRAY<T,VECTOR<int,3> >& u,ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat) const
{
    if(plan_u_to_u_hat_counts!=grid.Counts()){
        plan_u_to_u_hat_counts=grid.Counts();
        FFTW<T,3>::Destroy_Plan(plan_u_to_u_hat);
        plan_u_to_u_hat=FFTW<T,3>::Plan_R2C(grid.Counts(),u,u_hat);}
    FFTW<T,3>::Execute_R2C(plan_u_to_u_hat,u,u_hat);
}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class T> void FFT_3D<T>::
Inverse_Transform(ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat,ARRAY<T,VECTOR<int,3> >& u,bool normalize,bool preserve_u_hat) const
{
    if(plan_u_hat_to_u_counts!=grid.Counts()){
        plan_u_hat_to_u_counts=grid.Counts();
        FFTW<T,3>::Destroy_Plan(plan_u_hat_to_u);
        plan_u_hat_to_u=FFTW<T,3>::Plan_C2R(grid.Counts(),u_hat,u);}
    if(preserve_u_hat){
        u_hat_copy=u_hat;
        FFTW<T,3>::Execute_C2R(plan_u_hat_to_u,u_hat_copy,u);}
    else FFTW<T,3>::Execute_C2R(plan_u_hat_to_u,u_hat,u);
    if(normalize) u*=(T)1/grid.Counts().Product();
}
//#####################################################################

#else // NR-specific functions **********************************************************************************************************************************************************************************
#include <PhysBAM_Tools/Fourier_Transforms/FFT.h>
//#####################################################################
// Constructor
//#####################################################################
template<class T> FFT_3D<T>::
FFT_3D(const GRID<TV>& grid_input)
    :grid(grid_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FFT_3D<T>::
~FFT_3D()
{}
//#####################################################################
// Function Transform
//#####################################################################
template<class T> void FFT_3D<T>::
Transform(const ARRAY<T,VECTOR<int,3> >& u,ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat) const
{
    ARRAY<int> dim(grid.counts);
    data.Resize(2*grid.counts.x*grid.counts.y*grid.counts.z,false,false);
    int k=0;for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int ij=1;ij<=grid.counts.z;ij++){data(++k)=(float)u(i,j,ij);data(++k)=0;}
    NR_fourn(-1,dim,data);
    k=0;for(int i=0;i<=grid.counts.x-1;i++) for(int j=0;j<=grid.counts.y-1;j++) {for(int ij=0;ij<=grid.counts.z/2;ij++){u_hat(i,j,ij).re=data(++k);u_hat(i,j,ij).im=data(++k);} k+=grid.counts.z-2;}
}
//#####################################################################
// Function Inverse_Transform
//#####################################################################
template<class T> void FFT_3D<T>::
Inverse_Transform(ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat,ARRAY<T,VECTOR<int,3> >& u,bool normalize,bool preserve_u_hat) const
{
    ARRAY<int> dim(grid.counts);
    data.Resize(2*grid.counts.x*grid.counts.y*grid.counts.z,false,false);
    int k=0;
    for(int i=0;i<=grid.counts.x-1;i++){int negi=(i==0)?0:grid.counts.x-i;
        for(int j=0;j<=grid.counts.y-1;j++){int negj=(j==0)?0:grid.counts.y-j;
            for(int ij=0;ij<=grid.counts.z/2;ij++){data(++k)=(float)u_hat(i,j,ij).re;data(++k)=(float)u_hat(i,j,ij).im;}
            for(int ij=grid.counts.z/2+1;ij<=grid.counts.z-1;ij++){data(++k)=(float)u_hat(negi,negj,grid.counts.z-ij).re;data(++k)=-(float)u_hat(negi,negj,grid.counts.z-ij).im;}}}
    NR_fourn(+1,dim,data);
    if(normalize) {T coefficient=T(1.)/(grid.counts.x*grid.counts.y*grid.counts.z);k=0;for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int ij=1;ij<=grid.counts.z;ij++) {u(i,j,ij)=coefficient*data(++k);++k;}}
    else {k=0;for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int ij=1;ij<=grid.counts.z;ij++) {u(i,j,ij)=data(++k);++k;}}
}
//#####################################################################
#endif // ********************************************************************************************************************************************************************************************************

//#####################################################################
// Function Enforce_Real_Valued_Symmetry
//#####################################################################
// enforce symmetry so that the inverse transform is a real valued function
template<class T> void FFT_3D<T>::
Enforce_Real_Valued_Symmetry(ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat) const
{
    // imaginary part of the constant and cosine only terms are identically zero
    u_hat(0,0,0).im=u_hat(grid.counts.x/2,0,0).im=u_hat(0,grid.counts.y/2,0).im=u_hat(grid.counts.x/2,grid.counts.y/2,0).im=0;
    u_hat(0,0,grid.counts.z/2).im=u_hat(grid.counts.x/2,0,grid.counts.z/2).im=u_hat(0,grid.counts.y/2,grid.counts.z/2).im=u_hat(grid.counts.x/2,grid.counts.y/2,grid.counts.z/2).im=0;
    // enforce appropriate complex conjugates
    // front face - i.e. z=0
    for(int i=1;i<=grid.counts.x/2-1;i++){
        u_hat(grid.counts.x-i,0,0)=u_hat(i,0,0).Conjugated(); // reflect y=0 line
        u_hat(grid.counts.x-i,grid.counts.y/2,0)=u_hat(i,grid.counts.y/2,0).Conjugated();} // reflect y=grid.counts.y/2 line
    for(int j=1;j<=grid.counts.y/2-1;j++){u_hat(0,grid.counts.y-j,0)=u_hat(0,j,0).Conjugated();} // reflect x=0 line
    for(int i=1;i<=grid.counts.x-1;i++) for(int j=1;j<=grid.counts.y/2-1;j++){u_hat(grid.counts.x-i,grid.counts.y-j,0)=u_hat(i,j,0).Conjugated();} // reflect interior area
    // middle face - i.e. z=grid.counts.z/2
    for(int i=1;i<=grid.counts.x/2-1;i++){
        u_hat(grid.counts.x-i,0,grid.counts.z/2)=u_hat(i,0,grid.counts.z/2).Conjugated(); // reflect y=0 line
        u_hat(grid.counts.x-i,grid.counts.y/2,grid.counts.z/2)=u_hat(i,grid.counts.y/2,grid.counts.z/2).Conjugated();} // reflect y=grid.counts.y/2 line
    for(int j=1;j<=grid.counts.y/2-1;j++){u_hat(0,grid.counts.y-j,grid.counts.z/2)=u_hat(0,j,grid.counts.z/2).Conjugated();} // reflect x=0 line
    for(int i=1;i<=grid.counts.x-1;i++) for(int j=1;j<=grid.counts.y/2-1;j++){u_hat(grid.counts.x-i,grid.counts.y-j,grid.counts.z/2)=u_hat(i,j,grid.counts.z/2).Conjugated();} // reflect interior area
}
//#####################################################################
// Function First_Derivatives
//#####################################################################
template<class T> void FFT_3D<T>::
First_Derivatives(const ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat,ARRAY<COMPLEX<T> ,VECTOR<int,3> >& ux_hat,ARRAY<COMPLEX<T> ,VECTOR<int,3> >& uy_hat,ARRAY<COMPLEX<T> ,VECTOR<int,3> >& uz_hat) const
{
    VECTOR<T,3> coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    for(int i=0;i<=grid.counts.x-1;i++){T k1=coefficients.x*(i<=grid.counts.x/2?i:i-grid.counts.x);
        for(int j=0;j<=grid.counts.y-1;j++){T k2=coefficients.y*(j<=grid.counts.y/2?j:j-grid.counts.y);
            for(int ij=0;ij<=grid.counts.z/2;ij++){T k3=coefficients.z*ij;
                ux_hat(i,j,ij)=uy_hat(i,j,ij)=uz_hat(i,j,ij)=u_hat(i,j,ij).Rotated_Counter_Clockwise_90();
                ux_hat(i,j,ij)*=k1;uy_hat(i,j,ij)*=k2;uz_hat(i,j,ij)*=k3;}}}
    Enforce_Real_Valued_Symmetry(ux_hat);Enforce_Real_Valued_Symmetry(uy_hat);Enforce_Real_Valued_Symmetry(uz_hat);
}
//#####################################################################
// Make_Divergence_Free
//#####################################################################
// the constant fields are unchanged
template<class T> void FFT_3D<T>::
Make_Divergence_Free(ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat,ARRAY<COMPLEX<T> ,VECTOR<int,3> >& v_hat,ARRAY<COMPLEX<T> ,VECTOR<int,3> >& w_hat) const
{
    VECTOR<T,3> coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    // front face - i.e. z=0
    for(int i=1;i<=grid.counts.x/2;i++) u_hat(i,0,0)=COMPLEX<T>(0,0); // u=0 on the bottom
    for(int i=0;i<=grid.counts.x-1;i++){T k1=coefficients.x*(i<=grid.counts.x/2?i:i-grid.counts.x);
        for(int j=1;j<=grid.counts.y/2;j++){T k2=coefficients.y*j,one_over_k_squared=1/(sqr(k1)+sqr(k2));
            COMPLEX<T> correction=(k1*u_hat(i,j,0)+k2*v_hat(i,j,0))*one_over_k_squared;
            u_hat(i,j,0)-=correction*k1;v_hat(i,j,0)-=correction*k2;}}
    // volume
    for(int i=0;i<=grid.counts.x-1;i++){T k1=coefficients.x*(i<=grid.counts.x/2?i:i-grid.counts.x);
        for(int j=0;j<=grid.counts.y-1;j++){T k2=coefficients.y*(j<=grid.counts.y/2?j:j-grid.counts.y);
            for(int ij=1;ij<=grid.counts.z/2;ij++){T k3=coefficients.z*ij,one_over_k_squared=1/(sqr(k1)+sqr(k2)+sqr(k3));
                COMPLEX<T> correction=(k1*u_hat(i,j,ij)+k2*v_hat(i,j,ij)+k3*w_hat(i,j,ij))*one_over_k_squared;
                u_hat(i,j,ij)-=correction*k1;v_hat(i,j,ij)-=correction*k2;w_hat(i,j,ij)-=correction*k3;}}}
    Enforce_Real_Valued_Symmetry(u_hat);Enforce_Real_Valued_Symmetry(v_hat);Enforce_Real_Valued_Symmetry(w_hat);
}
//#####################################################################
// Function Filter_High_Frequencies
//#####################################################################
// Lanczos filter - doesn't change (0,0) frequency
template<class T> void FFT_3D<T>::
Filter_High_Frequencies(ARRAY<COMPLEX<T> ,VECTOR<int,3> >& u_hat,T scale) const
{
    T coefficient=2*T(pi)/sqrt((T)(sqr(grid.counts.x)+sqr(grid.counts.y)+sqr(grid.counts.z)));
    for(int i=1;i<=grid.counts.x/2;i++){ // skip the (0,0) case
        T temp=scale*coefficient*i,damping=sinc(temp);
        u_hat(i,0,0)*=damping;}
    for(int i=0;i<=grid.counts.x-1;i++){T i_frequency=T(i<=grid.counts.x/2 ? i:i-grid.counts.x);
        for(int j=1;j<=grid.counts.y/2;j++){T j_frequency=(T)j;
            T temp=scale*coefficient*sqrt(sqr(i_frequency)+sqr(j_frequency)),damping=sinc(temp);
            u_hat(i,j,0)*=damping;}}
    for(int i=0;i<=grid.counts.x-1;i++){T i_frequency=T(i<=grid.counts.x/2 ? i:i-grid.counts.x);
        for(int j=0;j<=grid.counts.y-1;j++){T j_frequency=T(j<=grid.counts.y/2 ? j:j-grid.counts.y);
            for(int ij=1;ij<=grid.counts.z/2;ij++){T ij_frequency=(T)ij;
                T temp=scale*coefficient*sqrt(sqr(i_frequency)+sqr(j_frequency)+sqr(ij_frequency)),damping=sinc(temp);
                u_hat(i,j,ij)*=damping;}}}
    Enforce_Real_Valued_Symmetry(u_hat);
}
//#####################################################################
template class FFT_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FFT_3D<double>;
#endif
