//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __GRID_ARRAYS_POLICY_DYADIC__
#define __GRID_ARRAYS_POLICY_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/ARRAYS_DYADIC_FORWARD.h>

namespace PhysBAM{

template<class T> class OCTREE_GRID;
template<class T> class QUADTREE_GRID;
template<class T> class BINTREE_GRID;

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T>
struct GRID_ARRAYS_POLICY<OCTREE_GRID<T> >
{
    typedef FLOOD_FILL_OCTREE<T> FLOOD_FILL;
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef ARRAY<T> FACE_ARRAYS;
    typedef ARRAY_BASE<T,ARRAY<T,int>,int> ARRAYS_BASE;
//#####################################################################
};

template<class T>
struct GRID_ARRAYS_POLICY<QUADTREE_GRID<T> >
{
    typedef FLOOD_FILL_QUADTREE<T> FLOOD_FILL;
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef ARRAY<T> FACE_ARRAYS;
    typedef ARRAY_BASE<T,ARRAY<T,int>,int> ARRAYS_BASE;
//#####################################################################
};

template<class T>
struct GRID_ARRAYS_POLICY<BINTREE_GRID<T> >
{
    typedef FLOOD_FILL_BINTREE<T> FLOOD_FILL;
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef ARRAY<T> FACE_ARRAYS;
    typedef ARRAY_BASE<T,ARRAY<T,int>,int> ARRAYS_BASE;
//#####################################################################
};

}
#endif
#endif
