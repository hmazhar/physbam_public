//#####################################################################
// Copyright 2009, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM__
#define __GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM__

namespace PhysBAM{
template<class T_GRID> struct COLLISION_GEOMETRY_COLLECTION_POLICY;
template<class T_GRID> struct GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV> class GRID;

template<class TV>
struct COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >
{
    typedef GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> > GRID_BASED_COLLISION_GEOMETRY;
};

}
#endif
