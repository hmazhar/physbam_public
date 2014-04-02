//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/AVERAGING_COLLIDABLE_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_RLE<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE:public ADVECTION<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:
    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    ARRAY<bool> &cell_valid_points_current,&cell_valid_points_next;
    T2 cell_crossover_replacement_value;
    bool extrapolate_to_revalidate_interpolation;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE<T_GRID,T2> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_RLE<T_GRID,T2,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_RLE<T_GRID,typename T_FACE_LOOKUP::NESTED_LOOKUP> velocity_averaging;
    AVERAGING_COLLIDABLE_RLE<T_GRID,T_FACE_LOOKUP> velocity_averaging_collidable;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,ARRAY<bool>& cell_valid_points_current_input,ARRAY<bool>& cell_valid_points_next_input,
        const T2& default_cell_replacement_value_input,const bool extrapolate_to_revalidate_interpolation_input)
        :body_list(body_list_input),cell_valid_points_current(cell_valid_points_current_input),cell_valid_points_next(cell_valid_points_next_input),
        cell_crossover_replacement_value(default_cell_replacement_value_input),extrapolate_to_revalidate_interpolation(extrapolate_to_revalidate_interpolation_input),
        linear_interpolation_collidable(body_list,&cell_valid_points_current,cell_crossover_replacement_value,extrapolate_to_revalidate_interpolation),
        velocity_averaging_collidable(body_list,0)
    {}

    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE()
    {}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2>& Z,const ARRAY<T2>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const ARRAY<T2>* Z_min_ghost,const ARRAY<T2>* Z_max_ghost,ARRAY<T2>* Z_min,ARRAY<T2>* Z_max)
    {assert(!Z_min && !Z_max);
    ARRAYS_COMPUTATIONS::Fill(cell_valid_points_next,true);
    for(CELL_ITERATOR iterator(grid,0);iterator;iterator++)if(iterator.Short()){int cell=iterator.Cell();
        T_BLOCK cell_block(grid,iterator.I()); 
        if(!cell_block || !body_list.Swept_Occupied_Block(cell_block)){
            cell_valid_points_next(cell)=true;
            TV sample_location=iterator.X()-dt*velocity_averaging.Face_To_Cell_Vector(iterator,face_velocities.Nested());
            T_BLOCK block(grid,sample_location);if(!block) continue;
            Z(cell)=linear_interpolation.From_Block_Cell(block,Z_ghost,sample_location);}
        else{
            TV grid_point_location=iterator.X(),length_and_direction=-dt*velocity_averaging_collidable.Face_To_Cell_Vector(cell_block,iterator,face_velocities),
                interpolation_point=grid_point_location+length_and_direction;
            if(body_list.Latest_Cell_Crossover(iterator,dt)){cell_valid_points_next(cell)=false;Z(cell)=linear_interpolation_collidable.default_cell_replacement_value;}
            else{
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(grid_point_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id))
                    interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                T_BLOCK block(grid,interpolation_point);if(!block){cell_valid_points_next(cell)=true;continue;}
                Z(cell)=linear_interpolation_collidable.From_Block_Cell(block,Z_ghost,interpolation_point,cell_valid_points_next(cell));}}}
    ARRAY<bool>::Exchange_Arrays(cell_valid_points_current,cell_valid_points_next);}

    void Average_To_Invalidated_Cells(const T_GRID& grid,const T2 default_value,ARRAY<T2>& values)
    {// average values collision aware in Gauss-Jacobi fashion
    bool done=false;ARRAY<PAIR<int,bool> > invalid_indices; // index and bool true if entry has been validated on iteration
    for(CELL_ITERATOR cell(grid,0);cell;cell++)if(cell.Short() && !cell_valid_points_current(cell.Cell())) invalid_indices.Append(PAIR<int,bool>(cell.Cell(),false));
    grid.Put_Ghost_Cells(false,cell_valid_points_current);
    const ARRAY<VECTOR<bool,T_GRID::dimension> >& cell_neighbors_visible=body_list.cell_neighbors_visible;
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& cell_neighbors=grid.Short_Cell_Neighbors();
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            int cell=invalid_indices(k).x;
            T2 sum=T2();int count=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                int neighbor1=cell_neighbors(cell)(2*axis-1),neighbor2=cell_neighbors(cell)(2*axis);
                if(neighbor1 && cell_neighbors_visible(neighbor1)(axis) && cell_valid_points_current(neighbor1)){sum+=values(neighbor1);count++;}
                if(neighbor2 && cell_neighbors_visible(cell)(axis) && cell_valid_points_current(neighbor2)){sum+=values(neighbor2);count++;}}
            if(count){values(cell)=sum/(T)count;invalid_indices(k).y=true;done=false;}}
        if(!done) for(int k=invalid_indices.m;k>=1;k--)if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
    // keep a copy of currently valid nodes (used for phi so we can revalidate the remaining nodes again after collision aware fast marching)
    // but important to initialize ghost nodes to true since currently cell_valid_points_current has them set to false
    cell_valid_points_next=cell_valid_points_current;
    grid.Put_Ghost_Cells(true,cell_valid_points_next);
    // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
    done=false;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            int cell=invalid_indices(k).x;
            T2 sum=T2();int count=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                int neighbor1=cell_neighbors(cell)(2*axis-1),neighbor2=cell_neighbors(cell)(2*axis);
                if(neighbor1){
                    if(cell_neighbors_visible(neighbor1)(axis)){if(cell_valid_points_current(neighbor1)){sum+=values(neighbor1);count++;}}
                    else{sum+=default_value;count++;}}
                if(neighbor2){
                    if(cell_neighbors_visible(cell)(axis)){if(cell_valid_points_current(neighbor2)){sum+=values(neighbor2);count++;}}
                    else{sum+=default_value;count++;}}}
            if(count){values(cell)=sum/(T)count;invalid_indices(k).y=true;done=false;}
            else values(cell)=default_value;}
        if(!done) for(int k=invalid_indices.m;k>=1;k--)if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
    grid.Put_Ghost_Cells(true,cell_valid_points_current);}

//#####################################################################
};
}
#endif
#endif
