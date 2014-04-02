//#####################################################################
// Copyright 2003-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PHILLIPS_SPECTRUM
//#####################################################################
#ifndef __PHILLIPS_SPECTRUM__
#define __PHILLIPS_SPECTRUM__

#include <PhysBAM_Tools/Fourier_Transforms/FFT_POLICY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class PHILLIPS_SPECTRUM
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_T;
    typedef typename T_ARRAYS_T::template REBIND<COMPLEX<T> >::TYPE T_ARRAYS_COMPLEX;
    typedef typename FFT_POLICY<TV>::FFT T_FFT;
public:
    GRID<TV> grid;
    int directional_factor_exponent;
    T length_cut_off,amplitude,gravity;
    TV wind;
    T_ARRAYS_COMPLEX h_initial;

    PHILLIPS_SPECTRUM(const GRID<TV>& grid_input)
        :grid(grid_input),directional_factor_exponent(6),length_cut_off((T)1),amplitude((T)1e-4),gravity((T)9.8),wind((T)10*TV::All_Ones_Vector())
    {}

//#####################################################################
    void Generate_Spectrum(RANDOM_NUMBERS<T>& random,bool verbose=false);
    void Get_H_Hat(T_ARRAYS_COMPLEX& h_hat);
    void Get_Separate_H_Hats(T_ARRAYS_COMPLEX& h_hat_1,T_ARRAYS_COMPLEX& h_hat_2);
    void Get_Heightfield(T_ARRAYS_T& h);
//#####################################################################
};
}
#endif
