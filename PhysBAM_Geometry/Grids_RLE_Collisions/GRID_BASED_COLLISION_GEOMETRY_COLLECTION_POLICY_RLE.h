//#####################################################################
// Copyright 2009, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_RLE__
#define __GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_RLE__

namespace PhysBAM{
template<class T_GRID> struct COLLISION_GEOMETRY_COLLECTION_POLICY;
template<class T_GRID> struct GRID_BASED_COLLISION_GEOMETRY_RLE;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;


template<class T>
struct COLLISION_GEOMETRY_COLLECTION_POLICY<RLE_GRID_2D<T> >
{
    typedef GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_2D<T> > GRID_BASED_COLLISION_GEOMETRY;
};

template<class T>
struct COLLISION_GEOMETRY_COLLECTION_POLICY<RLE_GRID_3D<T> >
{
    typedef GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_3D<T> > GRID_BASED_COLLISION_GEOMETRY;
};

}
#endif
#endif
