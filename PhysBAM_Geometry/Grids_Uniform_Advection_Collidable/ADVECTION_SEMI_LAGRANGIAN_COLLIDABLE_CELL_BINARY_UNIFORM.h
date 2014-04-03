#if 0
//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM:public ADVECTION<T_GRID,T2,T_FACE_LOOKUP>
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

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,T_ARRAYS_BOOL& cell_valid_points_current_input,
        T_ARRAYS_BOOL& cell_valid_points_next_input,const T2& default_cell_replacement_value_input,const bool extrapolate_to_revalidate_interpolation_input)
        :body_list(body_list_input),cell_valid_points_current(cell_valid_points_current_input),cell_valid_points_next(cell_valid_points_next_input),
        cell_crossover_replacement_value(default_cell_replacement_value_input),extrapolate_to_revalidate_interpolation(extrapolate_to_revalidate_interpolation_input),
        linear_interpolation_collidable(body_list,&cell_valid_points_current,cell_crossover_replacement_value,extrapolate_to_revalidate_interpolation),velocity_averaging_collidable(body_list,0)
    {}

    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM()
    {}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
    {PHYSBAM_NOT_IMPLEMENTED();
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();
        if(!body_list.Swept_Occupied_Cell_Center(cell)){
            TV interpolation_point=iterator.Location()-dt*velocity_averaging.Face_To_Cell_Vector(grid,cell,face_velocities.Nested());
            Z(cell)=linear_interpolation.Clamped_To_Array(grid,Z_ghost,interpolation_point);
            cell_valid_points_next(cell)=true;
            if(Z_min && Z_max){
                VECTOR<T2,2> extrema=linear_interpolation.Extrema_Clamped_To_Array(grid,*Z_min_ghost,*Z_max_ghost,interpolation_point);
                (*Z_min)(cell)=extrema.x;(*Z_max)(cell)=extrema.y;}}
        else{
            TV grid_point_location=iterator.Location(),length_and_direction=-dt*velocity_averaging_collidable.Face_To_Cell_Vector(grid,cell,face_velocities),
                interpolation_point=grid_point_location+length_and_direction;
            if(body_list.Latest_Cell_Crossover(cell,dt)){
                cell_valid_points_next(cell)=false;Z(cell)=linear_interpolation_collidable.default_cell_replacement_value;
                if(Z_min && Z_max) (*Z_min)(cell)=(*Z_max)(cell)=linear_interpolation_collidable.default_cell_replacement_value;}
            else{
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(grid_point_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id))
                    interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                Z(cell)=linear_interpolation_collidable.Clamped_To_Array(grid,Z_ghost,interpolation_point,cell_valid_points_next(cell));
                if(Z_min && Z_max){
                    VECTOR<T2,2> extrema=linear_interpolation_collidable.Extrema_Clamped_To_Array(grid,*Z_min_ghost,*Z_max_ghost,interpolation_point);
                    (*Z_min)(cell)=extrema.x;(*Z_max)(cell)=extrema.y;}}}}
    T_ARRAYS_BOOL::Exchange_Arrays(cell_valid_points_current,cell_valid_points_next);}

    T2 Compute_Revalidation_Value(const TV& from,const TV& to,const T2& current_invalid_value,const T2& default_value)
    {return default_value;}

    void Average_To_Invalidated_Cells(const T_GRID& grid,const T2 default_value,T_ARRAYS_T2& values)
    {// average values collision aware in Gauss-Jacobi fashion
    PHYSBAM_NOT_IMPLEMENTED();
    const T_ARRAYS_BOOL_DIMENSION& cell_neighbors_visible=body_list.cell_neighbors_visible;
    bool done=false;ARRAY<PAIR<TV_INT,bool> > invalid_indices; // index and bool true if entry has been validated on iteration
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(!cell_valid_points_current(iterator.Cell_Index()))
        invalid_indices.Append(PAIR<TV_INT,bool>(iterator.Cell_Index(),false));
    T_ARRAYS_BOOL::Put_Ghost(false,cell_valid_points_current,grid,3); // don't average from boundaries

    while(!done){done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            T2 sum=T2();int count=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                TV_INT min_cell=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_cell=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                if(cell_neighbors_visible(axis,min_cell) && cell_valid_points_current(min_cell)){sum+=values(min_cell);count++;}
                if(cell_neighbors_visible(axis,invalid_indices(k).x) && cell_valid_points_current(max_cell)){sum+=values(max_cell);count++;}}
            if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}

    // keep a copy of currently valid cells (used for phi so we can revalidate the remaining cells again after collision aware fast marching)
    // but important to initialize ghost cells to true since currently cell_valid_points_current has them set to false
    cell_valid_points_next=cell_valid_points_current;
    T_ARRAYS_BOOL::Put_Ghost(true,cell_valid_points_next,grid,3);

    // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
    done=false;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            T2 sum=T2();int count=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                TV_INT min_cell=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_cell=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                if(cell_neighbors_visible(axis,min_cell)){if(cell_valid_points_current(min_cell)){sum+=values(min_cell);count++;}}
                else{sum+=Compute_Revalidation_Value(grid.X(invalid_indices(k).x),grid.X(min_cell),values(invalid_indices(k).x),default_value);count++;}
                if(cell_neighbors_visible(axis,invalid_indices(k).x)){if(cell_valid_points_current(max_cell)){sum+=values(max_cell);count++;}}
                else{sum+=Compute_Revalidation_Value(grid.X(invalid_indices(k).x),grid.X(max_cell),values(invalid_indices(k).x),default_value);count++;}}
            if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}
            else values(invalid_indices(k).x)=default_value;}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
    T_ARRAYS_BOOL::Put_Ghost(true,cell_valid_points_current,grid,3);} // set valid for future advection

//#####################################################################
};
}
#endif
#endif
