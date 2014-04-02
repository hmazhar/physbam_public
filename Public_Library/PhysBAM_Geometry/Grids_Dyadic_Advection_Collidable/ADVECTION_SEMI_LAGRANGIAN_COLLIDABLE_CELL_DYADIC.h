//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/AVERAGING_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/AVERAGING_COLLIDABLE_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC:public ADVECTION<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT T_VECTOR_INT;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename T_GRID::CELL T_CELL;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;typedef typename T_GRID::UNIFORM_GRID T_UNIFORM_GRID;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:
    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    ARRAY<bool> &cell_valid_points_current,&cell_valid_points_next;
    T2 cell_crossover_replacement_value;
    bool extrapolate_to_revalidate_interpolation;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<T_GRID,T2> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_DYADIC<T_GRID,T2,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_DYADIC<T_GRID,typename T_FACE_LOOKUP::NESTED_LOOKUP> velocity_averaging;
    AVERAGING_COLLIDABLE_DYADIC<T_GRID,T_FACE_LOOKUP> velocity_averaging_collidable;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,ARRAY<bool>& cell_valid_points_current_input,ARRAY<bool>& cell_valid_points_next_input,
        const T2& default_cell_replacement_value_input,const bool extrapolate_to_revalidate_interpolation_input)
        :body_list(body_list_input),cell_valid_points_current(cell_valid_points_current_input),cell_valid_points_next(cell_valid_points_next_input),
        cell_crossover_replacement_value(default_cell_replacement_value_input),extrapolate_to_revalidate_interpolation(extrapolate_to_revalidate_interpolation_input),
        linear_interpolation_collidable(body_list,&cell_valid_points_current,cell_crossover_replacement_value,extrapolate_to_revalidate_interpolation),
        velocity_averaging_collidable(body_list,0)
    {}

    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC()
    {}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2>& Z,const ARRAY<T2>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const ARRAY<T2>* Z_min_ghost,const ARRAY<T2>* Z_max_ghost,ARRAY<T2>* Z_min,ARRAY<T2>* Z_max)
    {assert(!Z_min && !Z_max);
    ARRAY<T2> Z_node_ghost(grid.number_of_nodes,false);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Nodes(grid,Z_ghost,Z_node_ghost);
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,0);iterator.Valid();iterator.Next()){int cell=iterator.Cell_Index();
        if(!body_list.Swept_Occupied_Cell_Center(cell)){
            TV sample_location=iterator.Location()-dt*velocity_averaging.Face_To_Cell_Vector(grid,iterator.Cell_Pointer(),face_velocities.Nested());
            Z(cell)=linear_interpolation.From_Close_Cell_Cell(grid,iterator.Cell_Pointer(),Z_ghost,&Z_node_ghost,sample_location);
            cell_valid_points_next(cell)=true;}
        else{
            TV grid_point_location=iterator.Location(),length_and_direction=-dt*velocity_averaging_collidable.Face_To_Cell_Vector(grid,iterator.Cell_Pointer(),face_velocities),
                interpolation_point=grid_point_location+length_and_direction;
            if(body_list.Latest_Cell_Crossover(cell,dt)){cell_valid_points_next(cell)=false;Z(cell)=linear_interpolation_collidable.default_cell_replacement_value;}
            else{
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(grid_point_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id))
                    interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                T_CELL* base_cell=grid.Base_Cell_By_Neighbor_Path(iterator.Cell_Pointer(),interpolation_point);assert(base_cell);
                BLOCK_DYADIC<T_GRID> block(grid,base_cell);
                Z(cell)=linear_interpolation_collidable.From_Block_Cell(grid,block,Z_ghost,interpolation_point,cell_valid_points_next(cell));}}}
    ARRAY<bool>::Exchange_Arrays(cell_valid_points_current,cell_valid_points_next);}

private:
    void Get_Face_Values(const T_GRID& grid,T_CELL* cell,const ARRAY<bool>& cell_valid_points_current,const ARRAY<VECTOR<bool,T_GRID::dimension> >& cell_neighbors_visible,
        ARRAY<T_CELL*>& face_neighbors,const bool revalidate,ARRAY<T2>& values,T2& sum,T& weights,const T2& default_value)
    {ARRAY<VECTOR<T_CELL*,T_GRID::number_of_neighbors_per_cell> >& cell_neighbors=grid.Neighbors();
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        T_CELL* min_neighbor=cell_neighbors(cell->Cell())(2*axis-1),*max_neighbor=cell_neighbors(cell->Cell())(2*axis);
        if(min_neighbor){
            if(min_neighbor->Depth_Of_This_Cell()==grid.maximum_depth && cell_neighbors_visible(min_neighbor->Cell())(axis)){
                if(cell_valid_points_current(min_neighbor->Cell())){sum+=min_neighbor->Face_Size()*values(min_neighbor->Cell());weights+=min_neighbor->Face_Size();}
                else{sum+=min_neighbor->Face_Size()*default_value;weights+=min_neighbor->Face_Size();}}
            else{
                face_neighbors.Remove_All();
                cell->Get_All_Face_Neighbors(1,face_neighbors,&grid);
                for(int n=1;n<=face_neighbors.m;n++){
                    if(cell_valid_points_current(face_neighbors(n)->Cell())){sum+=face_neighbors(n)->Face_Size()*values(face_neighbors(n)->Cell());weights+=face_neighbors(n)->Face_Size();}
                    else if(revalidate){sum+=face_neighbors(n)->Face_Size()*default_value;weights+=face_neighbors(n)->Face_Size();}}}}
        if(max_neighbor){
            if(max_neighbor->Depth_Of_This_Cell()==grid.maximum_depth && cell_neighbors_visible(cell->Cell())(axis)){
                if(cell_valid_points_current(max_neighbor->Cell())){sum+=max_neighbor->Face_Size()*values(max_neighbor->Cell());weights+=max_neighbor->Face_Size();}
                else{sum+=max_neighbor->Face_Size()*default_value;weights+=max_neighbor->Face_Size();}}
            else{
                face_neighbors.Remove_All();
                cell->Get_All_Face_Neighbors(1,face_neighbors,&grid);
                for(int n=1;n<=face_neighbors.m;n++){
                    if(cell_valid_points_current(face_neighbors(n)->Cell())){sum+=face_neighbors(n)->Face_Size()*values(face_neighbors(n)->Cell());weights+=face_neighbors(n)->Face_Size();}
                    else if(revalidate){sum+=face_neighbors(n)->Face_Size()*default_value;weights+=face_neighbors(n)->Face_Size();}}}}}}
public:

    void Average_To_Invalidated_Cells(const T_GRID& grid,const T2 default_value,ARRAY<T2>& values)
    {// average values collision aware in Gauss-Jacobi fashion
    bool done=false;ARRAY<PAIR<T_CELL*,bool> > invalid_indices; // index and bool true if entry has been validated on iteration
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid);iterator.Valid();iterator.Next()) if(!cell_valid_points_current(iterator.Cell_Index()))
        invalid_indices.Append(PAIR<T_CELL*,bool>(iterator.Cell_Pointer(),false));
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,3,T_UNIFORM_GRID::GHOST_REGION);iterator.Valid();iterator.Next()) cell_valid_points_current(iterator.Cell_Index())=false;

    ARRAY<T_CELL*> face_neighbors;
    const ARRAY<VECTOR<bool,T_GRID::dimension> >& cell_neighbors_visible=body_list.cell_neighbors_visible;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            T_CELL* cell=invalid_indices(k).x;
            T2 sum=T2();T weights=0;
            Get_Face_Values(grid,invalid_indices(k).x,cell_valid_points_current,cell_neighbors_visible,face_neighbors,false,values,sum,weights,default_value);
            if(weights>0){values(cell->Cell())=sum/weights;invalid_indices(k).y=true;done=false;}}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x->Cell())=true;invalid_indices.Remove_Index_Lazy(k);}}
    // keep a copy of currently valid nodes (used for phi so we can revalidate the remaining nodes again after collision aware fast marching)
    // but important to initialize ghost nodes to true since currently cell_valid_points_current has them set to false
    cell_valid_points_next=cell_valid_points_current;
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,3,T_UNIFORM_GRID::GHOST_REGION);iterator.Valid();iterator.Next()) cell_valid_points_next(iterator.Cell_Index())=true;
    // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
    done=false;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            T_CELL* cell=invalid_indices(k).x;
            T2 sum=T2();T weights=0;
            Get_Face_Values(grid,invalid_indices(k).x,cell_valid_points_current,cell_neighbors_visible,face_neighbors,true,values,sum,weights,default_value);
            if(weights>0){values(cell->Cell())=sum/weights;invalid_indices(k).y=true;done=false;}
            else values(cell->Cell())=default_value;}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x->Cell())=true;invalid_indices.Remove_Index_Lazy(k);}}
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,3,T_UNIFORM_GRID::GHOST_REGION);iterator.Valid();iterator.Next()) cell_valid_points_current(iterator.Cell_Index())=true;}

//#####################################################################
};
}
#endif
#endif
