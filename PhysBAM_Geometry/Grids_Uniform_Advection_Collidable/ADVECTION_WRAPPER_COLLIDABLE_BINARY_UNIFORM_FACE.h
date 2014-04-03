//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_COLLIDABLE_BINARY_UNIFORM_FACE
//#####################################################################
#ifndef __ADVECTION_WRAPPER_COLLIDABLE_BINARY_UNIFORM_FACE__
#define __ADVECTION_WRAPPER_COLLIDABLE_BINARY_UNIFORM_FACE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_COLLIDABLE_FORWARD.h>
namespace PhysBAM{

// Assumes NESTED_ADVECTION is of type ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE
template<class T_GRID,class T2,class T_NESTED_LOOKUP,class T_NESTED_ADVECTION>
class ADVECTION_WRAPPER_COLLIDABLE_BINARY_UNIFORM_FACE:public ADVECTION<T_GRID,T2,T_NESTED_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_GRID::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename T_GRID::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<T_GRID,T_NESTED_LOOKUP> T_FACE_LOOKUP_COLLIDABLE;
    typedef typename T_GRID::BOUNDARY_SCALAR T_BOUNDARY;
public:
    T_NESTED_ADVECTION& nested_advection;
    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask;

    ADVECTION_WRAPPER_COLLIDABLE_BINARY_UNIFORM_FACE(T_NESTED_ADVECTION& nested_advection_input,T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
        :nested_advection(nested_advection_input),body_list(body_list_input),face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_NESTED_LOOKUP& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_NESTED_LOOKUP* Z_min_ghost,const T_NESTED_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {T_FACE_LOOKUP_COLLIDABLE Z_ghost_lookup(Z_ghost,body_list,&face_velocities_valid_mask);
    T_FACE_LOOKUP_COLLIDABLE V_lookup(face_velocities,body_list,&face_velocities_valid_mask);
    if(Z_min_ghost && Z_max_ghost){
        T_FACE_LOOKUP_COLLIDABLE Z_min_ghost_lookup(*Z_min_ghost,body_list,&face_velocities_valid_mask),Z_max_ghost_lookup(*Z_max_ghost,body_list,&face_velocities_valid_mask);
        nested_advection.Update_Advection_Equation_Face_Lookup(grid,Z,Z_ghost_lookup,V_lookup,boundary,dt,time,&Z_min_ghost_lookup,&Z_max_ghost_lookup,Z_min,Z_max);}
    else nested_advection.Update_Advection_Equation_Face_Lookup(grid,Z,Z_ghost_lookup,V_lookup,boundary,dt,time,0,0,Z_min,Z_max);}

//#####################################################################
};
}
#endif
