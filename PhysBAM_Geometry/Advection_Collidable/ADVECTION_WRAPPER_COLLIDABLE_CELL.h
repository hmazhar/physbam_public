//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_COLLIDABLE_CELL
//#####################################################################
#ifndef __ADVECTION_WRAPPER_COLLIDABLE_CELL__
#define __ADVECTION_WRAPPER_COLLIDABLE_CELL__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_DYADIC.h>
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_RLE.h>
#endif
namespace PhysBAM{

// Assumes NESTED_ADVECTION is of type ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_ 
template<class T_GRID,class T2,class T_NESTED_LOOKUP,class T_NESTED_ADVECTION,class T_FACE_LOOKUP_COLLIDABLE>
class ADVECTION_WRAPPER_COLLIDABLE_CELL:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename T_FACE_LOOKUP_COLLIDABLE::template REBIND_NESTED_LOOKUP<T_NESTED_LOOKUP>::TYPE T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:

    T_NESTED_ADVECTION& nested_advection;
    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask;

    ADVECTION_WRAPPER_COLLIDABLE_CELL(T_NESTED_ADVECTION& nested_advection_input,const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
        :nested_advection(nested_advection_input),body_list(body_list_input),face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
    {T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP V_lookup(face_velocities,body_list,&face_velocities_valid_mask);
    nested_advection.Update_Advection_Equation_Cell_Lookup(grid,Z,Z_ghost,V_lookup,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

//#####################################################################
};
}
#endif
