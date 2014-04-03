//#####################################################################
// Copyright 2005-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/RIGID_GEOMETRY_RASTERIZATION_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>::
GRID_BASED_COLLISION_GEOMETRY_DYADIC(T_GRID& grid_input)
    :GRID_BASED_COLLISION_GEOMETRY<T_GRID>(grid_input),use_collision_face_neighbors(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>::
~GRID_BASED_COLLISION_GEOMETRY_DYADIC()
{}
//##################################################################### 
// Function Initialize_Grids
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>::
Initialize_Grids()
{
    VECTOR<bool,T_GRID::dimension> all_true;all_true.Fill(true);
    cell_neighbors_visible.Resize(grid.number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(cell_neighbors_visible,all_true);
    face_neighbors_visible.Resize(grid.number_of_faces,false,false);ARRAYS_COMPUTATIONS::Fill(face_neighbors_visible,all_true);
}
//##################################################################### 
// Function Compute_Occupied_Blocks
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>::
Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor)
{
    ARRAY<bool>& occupied=with_body_motion?swept_occupied_blocks:occupied_blocks;
    occupied.Resize(grid.Block_Indices(grid.number_of_ghost_cells),false,false);ARRAYS_COMPUTATIONS::Fill(occupied,false);
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++)
        if(collision_geometry_collection.bodies(i) && collision_geometry_collection.bodies(i)->active)
            RASTERIZATION::Compute_Occupied_Blocks(*collision_geometry_collection.bodies(i),grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
}
//#####################################################################
// Function Compute_Grid_Visibility
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>::
Compute_Grid_Visibility()
{
    VECTOR<bool,T_GRID::dimension> all_true;all_true.Fill(true);
    cell_neighbors_visible.Resize(grid.number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(cell_neighbors_visible,all_true);
    face_neighbors_visible.Resize(grid.number_of_faces,false,false);ARRAYS_COMPUTATIONS::Fill(face_neighbors_visible,all_true);
    if(!collision_geometry_collection.bodies.m) return;

    // cell neighbors
    for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        T_CELL *cell1=iterator.First_Cell(),*cell2=iterator.Second_Cell();
        if(cell1->Depth_Of_This_Cell()==grid.maximum_depth&&cell2->Depth_Of_This_Cell()==grid.maximum_depth){
            ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1->Cell(),cell2->Cell(),collision_geometry_collection.bodies.m,objects);
            if(collision_geometry_collection.Intersection_Between_Points(cell1->Center(),cell2->Center(),&objects))cell_neighbors_visible(cell1->Cell())(iterator.Axis()+1)=false;}}

    // face neighbors
    for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        T_CELL* cell=iterator.Deepest_Cell();if(cell->Depth_Of_This_Cell()!=grid.maximum_depth)continue;
        for(int i=0;i<T_GRID::dimension;i++){
            int neighbor=iterator.Minimal_Neighbor(2*i+1);if(!neighbor)continue;
            DYADIC_GRID_ITERATOR_FACE<T_GRID> neighbor_iterator(grid,neighbor);
            ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell->Cell(),neighbor_iterator.Deepest_Cell_Index(),collision_geometry_collection.bodies.m,objects);
            if(collision_geometry_collection.Intersection_Between_Points(iterator.Location(),neighbor_iterator.Location(),&objects))face_neighbors_visible(iterator.Face_Index())(i+1)=false;}}
}
//#####################################################################
// Function Compute_Psi_N
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>::
Compute_Psi_N(ARRAY<bool>& psi_N,ARRAY<T>* face_velocities) const
{
    if(!collision_geometry_collection.bodies.m) return;
    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        T_CELL *cell1=iterator.Deepest_Cell(),*cell2=iterator.Other_Cell();int face_index=iterator.Face_Index();
        ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1->Cell(),cell2->Cell(),collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        COLLISION_GEOMETRY_ID body_id;int count=0;TV direction=cell2->Center()-cell1->Center();T length=(T).5*direction.Normalize(),velocity=0;
        RAY<TV> ray(cell1->Center(),direction,true);ray.t_max=length;ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity=collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[iterator.Axis()+1];}
        ray.Initialize(cell2->Center(),-direction,true);ray.t_max=length;ray.semi_infinite=false;
        if(collision_geometry_collection.Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity+=collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[iterator.Axis()+1];}
        if(count){psi_N(face_index)=true;if(face_velocities) (*face_velocities)(face_index)=velocity/count;}}
}
//#####################################################################
// Function Compute_Simplices_In_Cell
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>::
Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,INDEX>& simplices_in_cell,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies,
    int ghost_cells,T thickness,bool assume_active) const
{
    simplices_in_cell.Resize(grid.number_of_cells);
    for(COLLISION_GEOMETRY_ID i(1);i<=bodies.m;i++) if(assume_active || Is_Active(i)){ // Is_Active is not safe to call if different bodies list is used
        COLLISION_GEOMETRY<TV>* body=bodies(i);
        int n=body->Number_Of_Simplices();
        for(int e=1;e<=n;e++){
            RANGE<TV> box=body->World_Space_Simplex(e).Bounding_Box().Thickened(thickness);
            RANGE<TV_INT> bounding_grid_cells(grid.uniform_grid.Clamp_To_Cell(box.min_corner,ghost_cells+1),grid.uniform_grid.Clamp_To_Cell(box.max_corner,ghost_cells+1));
            for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next())
                simplices_in_cell(iterator.Cell_Index()).Append(PAIR<COLLISION_GEOMETRY_ID,int>(i,e));}}
}
//#####################################################################
#if COMPILE_WITH_BINTREE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY_DYADIC<BINTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY_DYADIC<BINTREE_GRID<double> >;
#endif
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<float> >;
template class GRID_BASED_COLLISION_GEOMETRY_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<double> >;
template class GRID_BASED_COLLISION_GEOMETRY_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif

#endif
