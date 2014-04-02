//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEEP_WATER_EVOLUTION
//#####################################################################
#ifndef __DEEP_WATER_EVOLUTION__
#define __DEEP_WATER_EVOLUTION__

#include <PhysBAM_Tools/Fourier_Transforms/FFT_1D.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFT_2D.h>
#include <PhysBAM_Tools/Fourier_Transforms/FFT_POLICY.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/PHILLIPS_SPECTRUM.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class DEEP_WATER_EVOLUTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,TV::dimension+1> TV_FULL;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_T;
    typedef typename T_ARRAYS_T::template REBIND<TV>::TYPE T_ARRAYS_TV;
    typedef typename T_ARRAYS_T::template REBIND<COMPLEX<T> >::TYPE T_ARRAYS_COMPLEX;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;
    typedef typename FFT_POLICY<TV>::FFT T_FFT;
public:
    T gravity;
    T period; // 0 for not periodic
    T lambda; // scaling for displacement (choppy waves)
    bool filter_high_frequencies;
    T high_frequency_cutoff;

    GRID<TV> grid;
    T_ARRAYS_T h; // height
    T_ARRAYS_TV Xh; // optional horizontal displacement
    T_ARRAYS_COMPLEX h_hat,h_hat1,h_hat2; // Fourier coefficients
    VECTOR<T_ARRAYS_COMPLEX,TV::m> dXh_hat; // Fourier coefficients
    VECTOR<T_ARRAYS_T,TV::m> displacement; // Real valued signals
    RANDOM_NUMBERS<T> random;
    PHILLIPS_SPECTRUM<TV> phillips_spectrum;
    T_FFT fft;

    bool use_surface_pressure;
    T_ARRAYS_T surface_pressure; // actually pressure/density (p = gh, not rho g h)

    TV_INT texture_cutoffs;

    DEEP_WATER_EVOLUTION()
        :gravity((T)9.8),period(0),lambda(0),filter_high_frequencies(false),high_frequency_cutoff((T).5),phillips_spectrum(grid),fft(grid),use_surface_pressure(false)
    {}

    void Use_Surface_Pressure(const bool use_surface_pressure_input=true)
    {use_surface_pressure=use_surface_pressure_input;if(use_surface_pressure) surface_pressure.Resize(grid.Domain_Indices());else surface_pressure.Clean_Memory();}

    void Filter_High_Frequencies(const T high_frequency_cutoff_input)
    {filter_high_frequencies=true;high_frequency_cutoff=high_frequency_cutoff_input;}

protected:
    T Omega(const T k) const
    {T omega=sqrt(gravity*k);
    if(period){
        T omega0=T(2*pi)/period;
        omega=((int)(omega/omega0))*omega0;}
    return omega;}
public:

//#####################################################################
    void Initialize();
    void Advance_Height(const T dt);
    void Get_Vertical_Velocity(T_ARRAYS_T& v) const;
    void Texture_Shallow_Water(const int frame,const STREAM_TYPE stream_type,const std::string& shallow_water_directory);
protected:
    void Set_H_Hats_From_Height();
    void Advance_H_Hats(const T dt);
//#####################################################################
};
}
#endif
