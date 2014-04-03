//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADVECTION_COLLIDABLE_POLICY_UNIFORM__
#define __ADVECTION_COLLIDABLE_POLICY_UNIFORM__

#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY_UNIFORM.h>

namespace PhysBAM{

template<class T_GRID>
struct ADVECTION_COLLIDABLE_POLICY
{
private:
    typedef typename T_GRID::SCALAR T;
public:
    // collidable
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM<T_GRID,T> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM<T_GRID> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE;
    // slip collidable
    //typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_SLIP_UNIFORM<T_GRID,T,typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE_SLIP> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_CELL;
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<T_GRID,typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE_SLIP> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE;
};

}
#endif
