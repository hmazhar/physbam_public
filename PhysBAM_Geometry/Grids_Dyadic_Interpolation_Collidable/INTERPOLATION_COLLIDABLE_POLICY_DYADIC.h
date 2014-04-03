//#####################################################################
// Copyright 2009, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __INTERPOLATION_COLLIDABLE_POLICY_DYADIC__
#define __INTERPOLATION_COLLIDABLE_POLICY_DYADIC__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_DYADIC_FORWARD.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
namespace PhysBAM{

template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;

template<class T>
struct INTERPOLATION_COLLIDABLE_POLICY<OCTREE_GRID<T> >
{
private:
    typedef OCTREE_GRID<T> T_GRID;
public:
    typedef AVERAGING_DYADIC<T_GRID> AVERAGING;
    // collidable
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID> FACE_LOOKUP_COLLIDABLE;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    // slip collidable
    typedef ARRAY<T> FACE_ARRAYS_SLIP;
    typedef FACE_LOOKUP_DYADIC<T_GRID> FACE_LOOKUP_SLIP;
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID> FACE_LOOKUP_COLLIDABLE_SLIP;
};

template<class T>
struct INTERPOLATION_COLLIDABLE_POLICY<QUADTREE_GRID<T> >
{
private:
    typedef QUADTREE_GRID<T> T_GRID;
public:
    typedef AVERAGING_DYADIC<T_GRID> AVERAGING;
    // collidable
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID> FACE_LOOKUP_COLLIDABLE;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    // slip collidable
    typedef ARRAY<T> FACE_ARRAYS_SLIP;
    typedef FACE_LOOKUP_DYADIC<T_GRID> FACE_LOOKUP_SLIP;
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID> FACE_LOOKUP_COLLIDABLE_SLIP;
};
//#####################################################################
}
#endif
#endif
