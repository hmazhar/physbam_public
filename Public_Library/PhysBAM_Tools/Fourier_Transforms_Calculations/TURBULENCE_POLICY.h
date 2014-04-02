//#####################################################################
// Copyright 2006, Geoffrey Irving, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TURBULENCE_POLICY 
//#####################################################################
#ifndef __TURBULENCE_POLICY__
#define __TURBULENCE_POLICY__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T> class TURBULENCE;
template<class T> class TURBULENCE_1D;
template<class T> class TURBULENCE_2D;
template<class T> class TURBULENCE_3D;

template<class TV> struct TURBULENCE_POLICY;

//#####################################################################
// 1D
//#####################################################################
template<class T>
struct TURBULENCE_POLICY<VECTOR<T,1> >
{
    typedef TURBULENCE_1D<T> TURBULENCE;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct TURBULENCE_POLICY<VECTOR<T,2> >
{
    typedef TURBULENCE_2D<T> TURBULENCE;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct TURBULENCE_POLICY<VECTOR<T,3> >
{
    typedef TURBULENCE_3D<T> TURBULENCE;
};
//#####################################################################
}
#endif
