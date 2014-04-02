//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM:public ADVECTION<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_TV;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:
    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    T_ARRAYS_BOOL &cell_valid_points_current,&cell_valid_points_next;
    T2 cell_crossover_replacement_value;
    bool extrapolate_to_revalidate_interpolation;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_UNIFORM<T_GRID,typename T_FACE_LOOKUP::NESTED_LOOKUP> velocity_averaging;
    AVERAGING_COLLIDABLE_UNIFORM<T_GRID,T_FACE_LOOKUP> velocity_averaging_collidable;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,T_ARRAYS_BOOL& cell_valid_points_current_input,
        T_ARRAYS_BOOL& cell_valid_points_next_input,const T2& default_cell_replacement_value_input,const bool extrapolate_to_revalidate_interpolation_input);
    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM();

    T2 Compute_Revalidation_Value(const TV& from,const TV& to,const T2& current_invalid_value,const T2& default_value)
    {return default_value;}

//#####################################################################
    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max);
    void Average_To_Invalidated_Cells(const T_GRID& grid,const T2 default_value,T_ARRAYS_T2& values);
//#####################################################################
};
}
#endif
