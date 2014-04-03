//#####################################################################
// Copyright 2003-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PHILLIPS_SPECTURM
//##################################################################### 
#include <PhysBAM_Tools/Fourier_Transforms/FFT_1D.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFT_2D.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/FREQUENCY_ITERATOR.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/PHILLIPS_SPECTRUM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Generate_Spectrum
//#####################################################################
template<class TV> void PHILLIPS_SPECTRUM<TV>::
Generate_Spectrum(RANDOM_NUMBERS<T>& random,bool verbose)
{
    TV wind_direction=wind.Normalized();
    T L=wind.Magnitude_Squared()/gravity;

    if(verbose){
        LOG::cout<<"Generating Phillips Spectrum"<<std::endl;
        LOG::cout<<"L = "<<L<<std::endl;
        LOG::cout<<"DX = "<<grid.dX<<std::endl;}

    TV coefficients=2*(T)pi/grid.domain.Edge_Lengths();

    typename FREQUENCY_ITERATOR<TV>::FULL_DOMAIN full_domain_type;
    h_initial.Resize(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid,full_domain_type));
    for(FREQUENCY_ITERATOR<TV> iterator(grid,full_domain_type);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
        TV k_vector=iterator.Frequency();
        T k=k_vector.Magnitude();

        T energy=0;
        if(k != 0){
            energy=amplitude*exp(-1/sqr(k*L))/pow(k,4+directional_factor_exponent); // power law
            // added abs in case directional_factor_exponent is odd
            energy*=abs(pow(TV::Dot_Product(k_vector,wind_direction),directional_factor_exponent)); // damps waves perpendicular to wind - does nothing in 1d but work as an amplitude
            energy*=exp(-sqr(k*length_cut_off));} // damps high frequency waves

        // find the fourier coeficients 
        T angle=random.Get_Uniform_Number((T)0,(T)2*(T)pi),r=random.Get_Gaussian();
        h_initial(I)=COMPLEX<T>::Polar(r*sqrt(energy),angle);}
}
//#####################################################################
// Function Get_H_Hat
//#####################################################################
template<class TV> void PHILLIPS_SPECTRUM<TV>::
Get_H_Hat(T_ARRAYS_COMPLEX& h_hat)
{
    for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
        h_hat(I)=h_initial(I)+h_initial(grid.counts-I).Conjugated();}
}
//#####################################################################
// Function Get_Separate_H_Hats
//#####################################################################
template<class TV> void PHILLIPS_SPECTRUM<TV>::
Get_Separate_H_Hats(T_ARRAYS_COMPLEX& h_hat1,T_ARRAYS_COMPLEX& h_hat2)
{
    PHYSBAM_ASSERT(T_ARRAYS_COMPLEX::Equal_Dimensions(h_hat1,h_hat2));
    for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
        h_hat1(I)=h_initial(I);
        h_hat2(I)=h_initial(grid.counts-I).Conjugated();}
}
//#####################################################################
// Function Get_Heightfield
//#####################################################################
template<class TV> void PHILLIPS_SPECTRUM<TV>::
Get_Heightfield(T_ARRAYS_T& h)
{
    PHYSBAM_ASSERT(h.Domain_Indices() == grid.Domain_Indices());
    T_ARRAYS_COMPLEX h_hat(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid));
    Get_H_Hat(h_hat);
    T_FFT fft(grid);
    fft.Inverse_Transform(h_hat,h,false,false);
}
//#####################################################################
template class PHILLIPS_SPECTRUM<VECTOR<float,1> >;
template class PHILLIPS_SPECTRUM<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PHILLIPS_SPECTRUM<VECTOR<double,1> >;
template class PHILLIPS_SPECTRUM<VECTOR<double,2> >;
#endif
