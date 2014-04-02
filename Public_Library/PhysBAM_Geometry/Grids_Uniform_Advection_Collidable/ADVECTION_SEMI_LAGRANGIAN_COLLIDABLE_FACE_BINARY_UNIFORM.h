//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM  
//##################################################################### 
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/AVERAGING_COLLIDABLE_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM:public ADVECTION<T_GRID,typename T_GRID::SCALAR,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;
    typedef ARRAY<bool,SIDED_FACE_INDEX<TV::dimension> > T_FACE_ARRAYS_SLIP_BOOL;
public:
    T_GRID_BASED_COLLISION_GEOMETRY& body_list;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM<T_GRID,T,T_FACE_LOOKUP> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_UNIFORM<T_GRID,typename T_FACE_LOOKUP::NESTED_LOOKUP> averaging;
    AVERAGING_COLLIDABLE_BINARY_UNIFORM<T_GRID,T_FACE_LOOKUP> averaging_collidable;
    T_FACE_ARRAYS_SLIP_BOOL& face_velocities_valid_mask;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM(T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input);
    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM();

    T Compute_Revalidation_Value(const int axis,const TV& from,const TV& to,const T& current_invalid_value,const T& default_value)
    {TV point;COLLISION_GEOMETRY_ID body_id;int aggregate_id;
    if(body_list.collision_geometry_collection.Intersection_Between_Points(from,to,body_id,aggregate_id,point)) return body_list.Object_Velocity(body_id,aggregate_id,point)[axis];
    else return default_value;}

//#####################################################################
    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max);
    void Average_To_Invalidated_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_values,T_FACE_ARRAYS_BOOL* faces_not_to_revalidate=0);
//#####################################################################
};
}
#endif
