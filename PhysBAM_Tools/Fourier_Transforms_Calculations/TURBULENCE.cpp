//#####################################################################
// Copyright 2002-2003, Robert Bridson, Ron Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms/FFT_2D.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFT_3D.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
using namespace PhysBAM;
//#####################################################################
// Function Generate_Random_Turbulence
//#####################################################################
template<class T> void TURBULENCE<T>::
Generate_Random_Turbulence(const GRID<VECTOR<T,2> >& grid,ARRAY<T,VECTOR<int,2> >& u,ARRAY<T,VECTOR<int,2> >& v) const
{
    int m=grid.counts.x,n=grid.counts.y;

    FFT_2D<T> fft(grid);
    ARRAY<COMPLEX<T> ,VECTOR<int,2> > u_hat(0,m-1,0,n/2),v_hat(0,m-1,0,n/2);
    
    T two_pi=(T)(2*pi);VECTOR<T,2> coefficients=two_pi/grid.domain.Edge_Lengths();
    for(int i=0;i<=m-1;i++) for(int j=0;j<=n/2;j++){
        T k1=coefficients.x*(i<=m/2 ? i:i-m),k2=coefficients.y*j,k=sqrt(sqr(k1)+sqr(k2));
        T area=two_pi*k; // circumference in 2D
        T energy=0;
        if(k > k_inertial) energy=constant*pow(epsilon,(T)two_thirds)*pow((T)k,-(T)five_thirds)/area;
        T r1=random->Get_Gaussian(),theta1=(T)pi*random->Get_Uniform_Number((T)0,(T)1); 
        T r2=random->Get_Gaussian(),theta2=(T)pi*random->Get_Uniform_Number((T)0,(T)1);
        T sqrt_energy_over_two=sqrt((T).5*energy);
        u_hat(i,j)=COMPLEX<T>::Polar(sqrt_energy_over_two*r1,theta1);
        v_hat(i,j)=COMPLEX<T>::Polar(sqrt_energy_over_two*r2,theta2);}
    
    fft.Enforce_Real_Valued_Symmetry(u_hat);fft.Enforce_Real_Valued_Symmetry(v_hat);
    if(incompressible) fft.Make_Divergence_Free(u_hat,v_hat);
    fft.Inverse_Transform(u_hat,u,false,false);fft.Inverse_Transform(v_hat,v,false,false);

    // rescale the final velocity
    if(rescaled_average_velocity){
        T average_velocity=0;
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) average_velocity+=sqrt(sqr(u(i,j))+sqr(v(i,j)));
        average_velocity/=(grid.counts.x*grid.counts.y);
        T scaling=rescaled_average_velocity/average_velocity;
        u*=scaling;v*=scaling;}
}
//#####################################################################
// Function Generate_Random_Turbulence
//#####################################################################
template<class T> void TURBULENCE<T>::
Generate_Random_Turbulence(const GRID<VECTOR<T,3> >& grid,ARRAY<T,VECTOR<int,3> >& u,ARRAY<T,VECTOR<int,3> >& v,ARRAY<T,VECTOR<int,3> >& w) const
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;

    FFT_3D<T> fft(grid);
    ARRAY<COMPLEX<T> ,VECTOR<int,3> > u_hat(0,m-1,0,n-1,0,mn/2),v_hat(0,m-1,0,n-1,0,mn/2),w_hat(0,m-1,0,n-1,0,mn/2);
    
    T four_pi=(T)(4*pi);VECTOR<T,3> coefficients=(T)(2*pi)/grid.domain.Edge_Lengths();
    for(int i=0;i<=m-1;i++) for(int j=0;j<=n-1;j++) for(int ij=0;ij<=mn/2;ij++){
        T k1=coefficients.x*(i<=m/2 ? i:i-m),k2=coefficients.y*(j<=n/2 ? j:j-n),k3=coefficients.z*ij;
        T k=sqrt(sqr(k1)+sqr(k2)+sqr(k3));
        T area=four_pi*sqr(k);
        T energy=0;
        if(k > k_inertial) energy=constant*pow(epsilon,(T)two_thirds)*pow((T)k,(T)-five_thirds)/area;
        T r1=random->Get_Gaussian(),theta1=(T)pi*random->Get_Uniform_Number((T)0,(T)1); 
        T r2=random->Get_Gaussian(),theta2=(T)pi*random->Get_Uniform_Number((T)0,(T)1);
        T r3=random->Get_Gaussian(),theta3=(T)pi*random->Get_Uniform_Number((T)0,(T)1);
        T sqrt_energy_over_two=sqrt((T).5*energy);
        u_hat(i,j,ij)=COMPLEX<T>::Polar(sqrt_energy_over_two*r1,theta1);
        v_hat(i,j,ij)=COMPLEX<T>::Polar(sqrt_energy_over_two*r2,theta2);
        w_hat(i,j,ij)=COMPLEX<T>::Polar(sqrt_energy_over_two*r3,theta3);}
    
    fft.Enforce_Real_Valued_Symmetry(u_hat);fft.Enforce_Real_Valued_Symmetry(v_hat);fft.Enforce_Real_Valued_Symmetry(w_hat);
    if(incompressible) fft.Make_Divergence_Free(u_hat,v_hat,w_hat);
    fft.Inverse_Transform(u_hat,u,false,false);fft.Inverse_Transform(v_hat,v,false,false);fft.Inverse_Transform(w_hat,w,false,false);

    // rescale the final velocity
    if(rescaled_average_velocity){
        T average_velocity=0;
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int ij=1;ij<=grid.counts.z;ij++)
            average_velocity+=sqrt(sqr(u(i,j,ij))+sqr(v(i,j,ij))+sqr(w(i,j,ij)));
        average_velocity/=(grid.counts.x*grid.counts.y*grid.counts.z);
        T scaling=rescaled_average_velocity/average_velocity;
        u*=scaling;v*=scaling;w*=scaling;}
}
//#####################################################################
template class TURBULENCE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TURBULENCE<double>;
#endif
