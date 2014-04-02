//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
using namespace PhysBAM;
//#####################################################################
// Function Interpolate_From_Nodes_To_Cells
//#####################################################################
template<class T_GRID> template<class T2> void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::
Interpolate_From_Nodes_To_Cells(const T_GRID& grid,const ARRAY<T2>& node_based,ARRAY<T2>& cell_based)
{
    ARRAY<T_CELL*>& cell_pointer_from_index=grid.Cell_Pointer_From_Index();cell_based.Resize(grid.number_of_cells);
    for(int i=1;i<=grid.number_of_cells;i++){T_CELL* cell=cell_pointer_from_index(i);if(cell){
        T2 sum=T2();for(int n=0;n<T_GRID::number_of_nodes_per_cell;n++)sum+=node_based(cell->Node(n));
        cell_based(i)=sum/(T)T_GRID::number_of_nodes_per_cell;}}
}
//#####################################################################
// Function Interpolate_From_Cell_Centers_To_Faces
//#####################################################################
template<class T_GRID> template<class T2> void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::
Interpolate_From_Cells_To_Faces(const T_GRID& grid,const ARRAY<T2>& cell_based,const ARRAY<T2>& node_based,ARRAY<T2>& face_based)
{
    for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
        T_CELL* deepest_cell=iterator.Deepest_Cell(),*other_cell=iterator.Other_Cell();
        if(other_cell&&deepest_cell->Depth_Of_This_Cell()==other_cell->Depth_Of_This_Cell())
            face_based(iterator.Face_Index())=(T).5*(cell_based(deepest_cell->Cell())+cell_based(other_cell->Cell()));
        else{
            face_based(iterator.Face_Index())=T2();
            for(int i=0;i<grid.number_of_nodes_per_face;i++)face_based(iterator.Face_Index())+=node_based(iterator.Face_Node(i));
            face_based(iterator.Face_Index())*=grid.one_over_number_of_nodes_per_face;}}
}
//#####################################################################
// Function Interpolate_From_Cells_To_Nodes
//#####################################################################
template<class T_GRID> template<class T2> void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::
Interpolate_From_Cells_To_Nodes(const T_GRID& grid,const ARRAY<T2>& cell_based,ARRAY<T2>& node_based)
{
    ARRAYS_COMPUTATIONS::Fill(node_based,T2());
    ARRAY<T> weight(grid.number_of_nodes);
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){T_CELL* cell=iterator.Cell_Pointer();
        T cell_weight=1/cell->DX().x;
        T2 weight_times_value=cell_weight*cell_based(iterator.Cell_Index());
        for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++){
            node_based(cell->Node(i))+=weight_times_value;
            weight(cell->Node(i))+=cell_weight;}}
    for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next())
        node_based(iterator.Node_Index())/=weight(iterator.Node_Index());
}
//#####################################################################
// Function Interpolate_From_Faces_To_Nodes
//#####################################################################
template<class T_GRID> void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::
Interpolate_From_Faces_To_Nodes(const T_GRID& grid,const ARRAY<T>& face_based,ARRAY<TV>& node_based)
{
    ARRAY<VECTOR<T,T_GRID::dimension> > weight(grid.number_of_nodes);
    ARRAYS_COMPUTATIONS::Fill(node_based,TV());
    for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
        T face_weight=1/iterator.Face_DX();
        T weight_times_value=face_weight*face_based(iterator.Face_Index());
        for(int i=0;i<grid.number_of_nodes_per_face;i++){
            node_based(iterator.Face_Node(i))[iterator.Axis()+1]+=weight_times_value;
            weight(iterator.Face_Node(i))(iterator.Axis()+1)+=face_weight;}}
    for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){int node=iterator.Node_Index();
        for(int axis=1;axis<=T_GRID::dimension;axis++)if(weight(node)(axis)) node_based(node)[axis]/=weight(node)(axis);}
}
//#####################################################################
// Function Interpolate_From_Faces_To_Nodes
//#####################################################################
template<class T_GRID> void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::
Interpolate_From_Faces_To_Nodes(const T_GRID& grid,const ARRAY<T>& face_based,ARRAY<T>& node_based)
{
    ARRAY<T> weight(grid.number_of_nodes);
    ARRAYS_COMPUTATIONS::Fill(node_based,T());
    for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
        T face_weight=1/iterator.Face_DX();
        T weight_times_value=face_weight*face_based(iterator.Face_Index());
        for(int i=0;i<grid.number_of_nodes_per_face;i++){
            node_based(iterator.Face_Node(i))+=weight_times_value;
            weight(iterator.Face_Node(i))+=face_weight;}}
    for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next())
        node_based(iterator.Node_Index())/=weight(iterator.Node_Index());
}
//#####################################################################
#define INSTANTIATION_HELPER(T,T_GRID) \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Nodes_To_Cells(const T_GRID& grid,const ARRAY<T_GRID::VECTOR_T>& node_based,ARRAY<T_GRID::VECTOR_T>& cell_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Nodes_To_Cells(const T_GRID& grid,const ARRAY<T>& node_based,ARRAY<T>& cell_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Cells_To_Faces(const T_GRID& grid,const ARRAY<T>& cell_based,const ARRAY<T>& node_based,ARRAY<T>& face_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Cells_To_Faces(const T_GRID& grid,const ARRAY<T_GRID::VECTOR_T>& cell_based,const ARRAY<T_GRID::VECTOR_T>& node_based,ARRAY<T_GRID::VECTOR_T>& face_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Cells_To_Nodes(const T_GRID& grid,const ARRAY<T>& cell_based,ARRAY<T>& node_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Cells_To_Nodes(const T_GRID& grid,const ARRAY<T_GRID::VECTOR_T>& cell_based,ARRAY<T_GRID::VECTOR_T>& node_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Cells_To_Nodes(const T_GRID& grid,const ARRAY<VECTOR<T_GRID::SCALAR,T_GRID::dimension+2> >& cell_based,ARRAY<VECTOR<T_GRID::SCALAR,T_GRID::dimension+2> >& node_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Faces_To_Nodes(const T_GRID& grid,const ARRAY<T>& face_based,ARRAY<TV>& node_based); \
    template void LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID >::Interpolate_From_Faces_To_Nodes(const T_GRID& grid,const ARRAY<T>& face_based,ARRAY<T>& node_based);
INSTANTIATION_HELPER(float,BINTREE_GRID<float>)
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
INSTANTIATION_HELPER(float,OCTREE_GRID<float>)
INSTANTIATION_HELPER(float,QUADTREE_GRID<float>)
#endif
#if !COMPILE_WITHOUT_DOUBLE_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
INSTANTIATION_HELPER(double,BINTREE_GRID<double>)
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
INSTANTIATION_HELPER(double,OCTREE_GRID<double>)
INSTANTIATION_HELPER(double,QUADTREE_GRID<double>)
#endif
#endif
#endif
