//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFT_POLICY 
//#####################################################################
#ifndef __FFT_POLICY__
#define __FFT_POLICY__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T> class FFT_1D;
template<class T> class FFT_2D;
template<class T> class FFT_3D;

template<class TV> struct FFT_POLICY;

//#####################################################################
// 1D
//#####################################################################
template<class T>
struct FFT_POLICY<VECTOR<T,1> >
{
    typedef FFT_1D<T> FFT;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct FFT_POLICY<VECTOR<T,2> >
{
    typedef FFT_2D<T> FFT;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct FFT_POLICY<VECTOR<T,3> >
{
    typedef FFT_3D<T> FFT;
};
//#####################################################################
}
#endif
