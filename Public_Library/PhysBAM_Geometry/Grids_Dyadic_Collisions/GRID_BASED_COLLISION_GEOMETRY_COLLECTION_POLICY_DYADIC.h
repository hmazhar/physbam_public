//#####################################################################
// Copyright 2009, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
// #if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_DYADIC__
#define __GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_DYADIC__

namespace PhysBAM{

template<class T_GRID> struct COLLISION_GEOMETRY_COLLECTION_POLICY;
template<class T_GRID> struct GRID_BASED_COLLISION_GEOMETRY_DYADIC;
template<class T> class BINTREE_GRID;
template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;

template<class T>
struct COLLISION_GEOMETRY_COLLECTION_POLICY<BINTREE_GRID<T> >
{
    typedef GRID_BASED_COLLISION_GEOMETRY_DYADIC<BINTREE_GRID<T> > GRID_BASED_COLLISION_GEOMETRY;
};

template<class T>
struct COLLISION_GEOMETRY_COLLECTION_POLICY<OCTREE_GRID<T> >
{
    typedef GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> > GRID_BASED_COLLISION_GEOMETRY;
};

template<class T>
struct COLLISION_GEOMETRY_COLLECTION_POLICY<QUADTREE_GRID<T> >
{
    typedef GRID_BASED_COLLISION_GEOMETRY_DYADIC<QUADTREE_GRID<T> > GRID_BASED_COLLISION_GEOMETRY;
};

}
#endif
#endif
