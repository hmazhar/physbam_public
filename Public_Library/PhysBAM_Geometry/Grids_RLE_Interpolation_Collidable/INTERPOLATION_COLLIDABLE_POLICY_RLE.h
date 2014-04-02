//#####################################################################
// Copyright 2009, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __INTERPOLATION_COLLIDABLE_POLICY_RLE__
#define __INTERPOLATION_COLLIDABLE_POLICY_RLE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_RLE_FORWARD.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
namespace PhysBAM{

template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;
template<class T,int d> class VECTOR;

template<class T>
struct INTERPOLATION_COLLIDABLE_POLICY<RLE_GRID_2D<T> >
{
private:
    typedef RLE_GRID_2D<T> T_GRID;
    typedef VECTOR<T,2> TV;
public:
    // normal
    typedef FACE_LOOKUP_RLE<T_GRID> FACE_LOOKUP;
    typedef AVERAGING_RLE<T_GRID> AVERAGING;
    typedef INTERPOLATION_RLE<T_GRID,T> INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T> LINEAR_INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,TV> LINEAR_INTERPOLATION_VECTOR;
    // collidable
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID,FACE_LOOKUP> FACE_LOOKUP_COLLIDABLE;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    // slip collidable
    typedef ARRAY<T> FACE_ARRAYS_SLIP;
    typedef FACE_LOOKUP_RLE<T_GRID> FACE_LOOKUP_SLIP;
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID> FACE_LOOKUP_COLLIDABLE_SLIP;
};

template<class T>
struct INTERPOLATION_COLLIDABLE_POLICY<RLE_GRID_3D<T> >
{
private:
    typedef RLE_GRID_3D<T> T_GRID;
    typedef VECTOR<T,3> TV;
public:
    // normal
    typedef FACE_LOOKUP_RLE<T_GRID> FACE_LOOKUP;
    typedef AVERAGING_RLE<T_GRID> AVERAGING;
    typedef INTERPOLATION_RLE<T_GRID,T> INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T> LINEAR_INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,TV> LINEAR_INTERPOLATION_VECTOR;
    // collidable
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID,FACE_LOOKUP> FACE_LOOKUP_COLLIDABLE;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    // slip collidable
    typedef ARRAY<T> FACE_ARRAYS_SLIP;
    typedef FACE_LOOKUP_RLE<T_GRID> FACE_LOOKUP_SLIP;
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID> FACE_LOOKUP_COLLIDABLE_SLIP;
};
//#####################################################################
}
#endif
#endif
