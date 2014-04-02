//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYADIC_GRID_POLICY
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __DYADIC_GRID_POLICY__
#define __DYADIC_GRID_POLICY__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
namespace PhysBAM{

template<class TV> struct DYADIC_TAG{};

template<class T> class BINTREE_GRID;
template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;

template<class TV> struct DYADIC_GRID_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct DYADIC_GRID_POLICY<VECTOR<T,0> >
{
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct DYADIC_GRID_POLICY<VECTOR<T,1> >
{
    typedef BINTREE_GRID<T> DYADIC_GRID;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct DYADIC_GRID_POLICY<VECTOR<T,2> >
{
    typedef QUADTREE_GRID<T> DYADIC_GRID;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct DYADIC_GRID_POLICY<VECTOR<T,3> >
{
    typedef OCTREE_GRID<T> DYADIC_GRID;
};
//#####################################################################

}
#endif
#endif
