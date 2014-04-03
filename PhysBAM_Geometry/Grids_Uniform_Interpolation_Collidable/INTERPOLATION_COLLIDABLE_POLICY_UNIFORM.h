//#####################################################################
// Copyright 2009, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERPOLATION_COLLIDABLE_POLICY_UNIFORM__
#define __INTERPOLATION_COLLIDABLE_POLICY_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
namespace PhysBAM{

template<class T_GRID>
struct INTERPOLATION_COLLIDABLE_POLICY
{
private:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;
public:
    typedef AVERAGING_UNIFORM<T_GRID> AVERAGING;
    typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> FACE_LOOKUP_COLLIDABLE;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM<T_GRID,T> LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    // slip collidable
    typedef ARRAY<T,SIDED_FACE_INDEX<TV::dimension> > FACE_ARRAYS_SLIP;
    typedef FACE_LOOKUP_UNIFORM<T_GRID> FACE_LOOKUP_SLIP;
    typedef FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID> FACE_LOOKUP_COLLIDABLE_SLIP;
};
//#####################################################################
}
#endif
