//#####################################################################
// Copyright 2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Computations/RIGID_GEOMETRY_RASTERIZATION_RLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
GRID_BASED_COLLISION_GEOMETRY_RLE(T_GRID& grid_input)
    :GRID_BASED_COLLISION_GEOMETRY<T_GRID>(grid_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
~GRID_BASED_COLLISION_GEOMETRY_RLE()
{}
//##################################################################### 
// Function Initialize_Grids
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
Initialize_Grids()
{
    cell_neighbors_visible.Resize(grid.number_of_cells,false,false);cell_neighbors_visible.Fill(true);
    face_neighbors_visible.Resize(grid.number_of_faces,false,false);face_neighbors_visible.Fill(true);
}
//##################################################################### 
// Function Compute_Occupied_Blocks
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor)
{
    ARRAY<bool>& occupied=with_body_motion?swept_occupied_blocks:occupied_blocks;
    occupied.Resize(grid.number_of_blocks,false,false);occupied.Fill(false);
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++)
        if(collision_geometry_collection.bodies(i) && collision_geometry_collection.bodies(i)->active)
            RASTERIZATION::Compute_Occupied_Blocks(*collision_geometry_collection.bodies(i),grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout<<(with_body_motion?"swept ":"")<<"occupied cells = "<<occupied.Number_True()<<" out of "<<occupied.m<<std::endl;
#endif
}
//#####################################################################
// Function Compute_Grid_Visibility
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
Compute_Grid_Visibility()
{
    cell_neighbors_visible.Resize(grid.number_of_cells,false,false);cell_neighbors_visible.Fill(true);
    face_neighbors_visible.Resize(grid.number_of_faces,false,false);face_neighbors_visible.Fill(true);
    if(!collision_geometry_collection.bodies.m) return;

    // cell neighbors
    T_GRID::template Face_Loop<Compute_Cell_Neighbors_Visible_Helper>(*this);

    // face neighbors
    ARRAY<int> face_cell(grid.number_of_faces);
    for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){if(cell.Short()) for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++)face_cell(cell.Face(i))=cell.Cell();}
    T_GRID::template Face_Loop<Compute_Face_Neighbors_Visible_Helper>(*this,face_cell);
}
//#####################################################################
// Function Compute_Cell_Neighbors_Visible_Helper
//#####################################################################
template<class T_GRID> template<class T_FACE> void GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
Compute_Cell_Neighbors_Visible_Helper::Apply(GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list)
{
    int axis=T_FACE::Axis();
    T_GRID& grid=body_list.grid;
    for(T_FACE face(grid,grid.number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){
        ARRAY<COLLISION_GEOMETRY_ID> objects;body_list.objects_in_cell.Get_Objects_For_Cells(face.cell1.Cell(),face.cell2.Cell(),body_list.collision_geometry_collection.bodies.m,objects);
        if(body_list.collision_geometry_collection.Intersection_Between_Points(face.cell1.X(),face.cell2.X(),&objects)) body_list.cell_neighbors_visible(face.cell1.Cell())(axis)=false;}
}
//#####################################################################
// Function Compute_Face_Neighbors_Visible_Helper
//#####################################################################
template<class T_GRID> template<class T_FACE> void GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
Compute_Face_Neighbors_Visible_Helper::Apply(GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list,ARRAY<int>& face_cell)
{
    T_GRID& grid=body_list.grid;
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& face_neighbors=grid.Short_Face_Neighbors();
    for(T_FACE face(grid,grid.number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){
        int f=face.Face();TV X=face.X();
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            int neighbor=face_neighbors(f)(2*axis);if(!neighbor) continue;
            ARRAY<COLLISION_GEOMETRY_ID> objects;body_list.objects_in_cell.Get_Objects_For_Cells(face_cell(f),face_cell(neighbor),body_list.collision_geometry_collection.bodies.m,objects);
            TV neighbor_X=X;neighbor_X[axis]+=grid.uniform_grid.dX[axis];
            if(body_list.collision_geometry_collection.Intersection_Between_Points(X,neighbor_X,&objects)) body_list.face_neighbors_visible(f)(axis)=false;}}
}
//#####################################################################
// Function Compute_Psi_N
//#####################################################################
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
Compute_Psi_N(ARRAY<bool>& psi_N,ARRAY<T>* face_velocities) const
{
    T_GRID::template Face_Loop<Compute_Psi_N_Helper>(*this,psi_N,face_velocities);
}
//##################################################################### 
// Function Compute_Psi_N_Helper
//##################################################################### 
template<class T_GRID> template<class T_FACE> void GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>::
Compute_Psi_N_Helper::Apply(const GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list,ARRAY<bool>& psi_N,ARRAY<T>* face_velocities)
{
    const int axis=T_FACE::Axis();
    const T_GRID& grid=body_list.grid;
    for(T_FACE face(grid,0);face;face++)if(face.Both_Cells_Short()){
        int f=face.Face(),c1=face.cell1.Cell(),c2=face.cell2.Cell();
        ARRAY<COLLISION_GEOMETRY_ID> objects;body_list.objects_in_cell.Get_Objects_For_Cells(c1,c2,body_list.collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
        COLLISION_GEOMETRY_ID body_id;int count=0;T velocity=0;
        RAY<TV> ray(face.cell1.X(),TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.uniform_grid.dX[axis];ray.semi_infinite=false;
        if(body_list.Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity=body_list.collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[axis];}
        ray.Initialize(face.cell2.X(),-TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.uniform_grid.dX[axis];ray.semi_infinite=false;
        if(body_list.Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)){count++;
            velocity+=body_list.collision_geometry_collection(body_id).Pointwise_Object_Pseudo_Velocity(ray.aggregate_id,ray.Point(ray.t_max),
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE)[axis];}
        if(count){psi_N(f)=true;if(face_velocities) (*face_velocities)(f)=velocity/count;}}
}
//#####################################################################
template class GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_2D<float> >;
template class GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_2D<double> >;
template class GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_3D<double> >;
#endif
#endif
