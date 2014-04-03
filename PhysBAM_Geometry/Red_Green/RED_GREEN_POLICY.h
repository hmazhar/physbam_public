//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_POLICY 
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __RED_GREEN_POLICY__
#define __RED_GREEN_POLICY__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
namespace PhysBAM{

template<class T> class RED_TRIANGLE;
template<class T> class RED_TETRAHEDRON;
template<class T> class RED_GREEN_GRID_2D;
template<class T> class RED_GREEN_GRID_3D;

template<class TV> struct RED_GREEN_POLICY;

//#####################################################################
// 2D
//#####################################################################
template<class T>
struct RED_GREEN_POLICY<VECTOR<T,2> >
{
    typedef RED_TRIANGLE<T> RED_SIMPLEX;
    typedef RED_GREEN_GRID_2D<T> GRID_T;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct RED_GREEN_POLICY<VECTOR<T,3> >
{
    typedef RED_TETRAHEDRON<T> RED_SIMPLEX;
    typedef RED_GREEN_GRID_3D<T> GRID_T;
};
//#####################################################################
}
#endif
#endif
