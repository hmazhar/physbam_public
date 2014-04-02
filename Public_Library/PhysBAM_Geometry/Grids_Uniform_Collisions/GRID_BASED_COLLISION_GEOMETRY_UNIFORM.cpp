//#####################################################################
// Copyright 2005-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_POINT_SIMPLEX_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_POINT_SIMPLEX_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/RIGID_GEOMETRY_RASTERIZATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
GRID_BASED_COLLISION_GEOMETRY_UNIFORM(T_GRID& grid_input)
    :GRID_BASED_COLLISION_GEOMETRY<T_GRID>(grid_input),use_collision_face_neighbors(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
~GRID_BASED_COLLISION_GEOMETRY_UNIFORM()
{}
//##################################################################### 
// Function Initialize_Grids
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
Initialize_Grids()
{
    collision_thickness=(T)1e-3*grid.Minimum_Edge_Length();
    collision_geometry_collection.Set_Collision_Body_Thickness(collision_thickness);

    assert(grid.Is_MAC_Grid());
    TV_BOOL all_true;all_true.Fill(true);
    cell_neighbors_visible.Resize(grid.Domain_Indices(3),false);cell_neighbors_visible.Fill(all_true); // initialize here so collision aware redistancing works in Initialize
    face_neighbors_visible.Resize(grid,1,false);face_neighbors_visible.Fill(all_true);
}
//##################################################################### 
// Function Compute_Occupied_Blocks
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor)
{
    T_ARRAYS_BOOL& occupied=with_body_motion?swept_occupied_blocks:occupied_blocks;
    occupied.Resize(grid.Block_Indices(3),false,false);occupied.Fill(false);
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++)
        if(collision_geometry_collection.bodies(i) && collision_geometry_collection.bodies(i)->active)
            RASTERIZATION::Compute_Occupied_Blocks(*collision_geometry_collection.bodies(i),grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
}
//#####################################################################
// Function Compute_Grid_Visibility
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
Compute_Grid_Visibility()
{
    TV_BOOL all_true;all_true.Fill(true);
    cell_neighbors_visible.Fill(all_true);face_neighbors_visible.Fill(all_true);
    if(!collision_geometry_collection.bodies.m) return;

    // cell neighbors
    for(FACE_ITERATOR iterator(grid,3,T_GRID::INTERIOR_REGION);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1,cell2,collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        if(collision_geometry_collection.Intersection_Between_Points(grid.Center(cell1),grid.Center(cell2),&objects)) cell_neighbors_visible(cell1)(iterator.Axis())=false;}

    // face neighbors
    for(FACE_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
        for(int direction=1;direction<=T_GRID::dimension;direction++){TV_INT direction_offset=TV_INT::Axis_Vector(direction);
            ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(iterator.Second_Cell_Index(),iterator.Second_Cell_Index()+direction_offset,collision_geometry_collection.bodies.m,objects);
            if(!objects.m) continue;
            if(collision_geometry_collection.Intersection_Between_Points(iterator.Location(),grid.Face(axis,iterator.Face_Index()+direction_offset),&objects))
                face_neighbors_visible.Component(axis)(iterator.Face_Index())(direction)=false;}}
}
//#####################################################################
// Function Compute_Psi_N
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
Compute_Psi_N(T_FACE_ARRAYS_BOOL& psi_N,T_FACE_ARRAYS_SCALAR* face_velocities) const
{
    if(!collision_geometry_collection.bodies.m) return;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1,cell2,collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();COLLISION_GEOMETRY_ID body_id;int count=0;T velocity=0;
        RAY<TV> ray(iterator.First_Cell_Center(),TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity=collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[axis];}
        ray.Initialize(iterator.Second_Cell_Center(),-TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity+=collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[axis];}
        if(count){psi_N.Component(axis)(face_index)=true;if(face_velocities) (*face_velocities).Component(axis)(face_index)=velocity/count;}}
}
//#####################################################################
// Function Compute_Psi_N
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
Compute_Psi_N_Zero_Velocity(T_FACE_ARRAYS_BOOL& psi_N,T_FACE_ARRAYS_SCALAR* face_velocities) const
{
    if(!collision_geometry_collection.bodies.m) return;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1,cell2,collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();COLLISION_GEOMETRY_ID body_id;int count=0;
        RAY<TV> ray(iterator.First_Cell_Center(),TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)) count++;
        ray.Initialize(iterator.Second_Cell_Center(),-TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
        if(Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)) count++;
        if(count){psi_N.Component(axis)(face_index)=true;if(face_velocities) (*face_velocities)(axis,face_index)=(T)0;}}
}
//#####################################################################
// Function Compute_Simplices_In_Cell
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>::
Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT>& simplices_in_cell,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies,
    int ghost_cells,T thickness,bool assume_active) const
{
    simplices_in_cell.Resize(grid.Domain_Indices(ghost_cells+1));
    for(COLLISION_GEOMETRY_ID i(1);i<=bodies.m;i++) if(assume_active || Is_Active(i)){ // Is_Active is not safe to call if different bodies list is used
        COLLISION_GEOMETRY<TV>* body=bodies(i);
        int n=body->Number_Of_Simplices();
        for(int e=1;e<=n;e++){
            BOX<TV> box=body->World_Space_Simplex(e).Bounding_Box().Thickened(thickness);
            RANGE<TV_INT> bounding_grid_cells(grid.Clamp_To_Cell(box.min_corner,ghost_cells+1),grid.Clamp_To_Cell(box.max_corner,ghost_cells+1));
            for(CELL_ITERATOR iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next())
                simplices_in_cell(iterator.Cell_Index()).Append(PAIR<COLLISION_GEOMETRY_ID,int>(i,e));}}
}
//#####################################################################
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<float,1> > >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<float,2> > >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<double,1> > >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<double,2> > >;
template class GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
