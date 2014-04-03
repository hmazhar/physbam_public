//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_RLE_MPI.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Data_Structures/GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#endif
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_RLE_MPI<T_GRID>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&)
{
    int m=local_grid.uniform_grid.counts.x;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,1> >(1,m-1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,1> >(0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,1> >(m,m),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    this->template Find_Boundary_Indices_In_Region<typename T_GRID::FACE_X_ITERATOR>(1,RANGE<VECTOR<int,1> >(1,1),cell_index_to_matrix_index);
    this->template Find_Boundary_Indices_In_Region<typename T_GRID::FACE_X_ITERATOR>(2,RANGE<VECTOR<int,1> >(m-1,m-1),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_RLE_MPI<T_GRID>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&)
{
    int m=local_grid.uniform_grid.counts.x,mn=local_grid.uniform_grid.counts.z;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,2> >(1,m-1,1,mn-1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,2> >(0,0,1,mn-1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,2> >(m,m,1,mn-1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(3,RANGE<VECTOR<int,2> >(1,m-1,0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(4,RANGE<VECTOR<int,2> >(1,m-1,mn,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    this->template Find_Boundary_Indices_In_Region<typename T_GRID::FACE_X_ITERATOR>(1,RANGE<VECTOR<int,2> >(1,1,1,mn-1),cell_index_to_matrix_index);
    this->template Find_Boundary_Indices_In_Region<typename T_GRID::FACE_X_ITERATOR>(2,RANGE<VECTOR<int,2> >(m-1,m-1,1,mn-1),cell_index_to_matrix_index);
    this->template Find_Boundary_Indices_In_Region<typename T_GRID::FACE_Z_ITERATOR>(3,RANGE<VECTOR<int,2> >(1,m-1,1,1),cell_index_to_matrix_index);
    this->template Find_Boundary_Indices_In_Region<typename T_GRID::FACE_Z_ITERATOR>(4,RANGE<VECTOR<int,2> >(1,m-1,mn-1,mn-1),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices_In_Region
//#####################################################################
template<class T_GRID> void LAPLACE_RLE_MPI<T_GRID>::
Find_Matrix_Indices_In_Region(const int region_index,const T_BOX_HORIZONTAL_INT& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,
    ARRAY<ARRAY<int> >& matrix_index_to_cell_index_array)
{
    if(region_index) for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).min_corner=filled_region_cell_count(color)+1;
    else for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).interior_indices.min_corner=filled_region_cell_count(color)+1;
    for(RLE_GRID_SIMPLE_ITERATOR<T_GRID,CELL_ITERATOR> cell(local_grid,region);cell;cell++){int c=cell.index;
        int color=filled_region_colors(c);if(color<1) continue;
        int new_index=++filled_region_cell_count(color);cell_index_to_matrix_index(c)=new_index;
        matrix_index_to_cell_index_array(color)(new_index)=c;}
    if(region_index) for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).max_corner=filled_region_cell_count(color);
    else for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).interior_indices.max_corner=filled_region_cell_count(color);
}
//#####################################################################
// Function Find_Boundary_Indices_In_Region
//#####################################################################
template<class T_GRID> template<class T_FACE> void LAPLACE_RLE_MPI<T_GRID>::
Find_Boundary_Indices_In_Region(const int side,const T_BOX_HORIZONTAL_INT& region,ARRAY<int>& cell_index_to_matrix_index)
{
    int horizontal_axis=T_FACE::Horizontal_Axis(),cell_side=side&1;assert(2*horizontal_axis-1<=side && side<=2*horizontal_axis);
    T_BOX_HORIZONTAL_INT face_region=region;if(!cell_side) face_region+=TV_HORIZONTAL_INT::Axis_Vector(horizontal_axis);
    // count boundary indices
    ARRAY<int> counts(partitions.m);
    int last_cell=0;
    for(T_FACE face(local_grid,face_region);face;face++)if(!psi_N(face.Face())){
        CELL_ITERATOR& cell=cell_side?face.cell2:face.cell1;
        int c=cell.Cell(),color=filled_region_colors(cell.Cell());
        if(c!=last_cell && counts.Valid_Index(color)){last_cell=c;counts(color)++;if(cell.Long()) counts(color)++;}}
    // fill boundary indices
    for(int color=1;color<=partitions.m;color++)partitions(color).boundary_indices(side).Resize(counts(color));
    ARRAYS_COMPUTATIONS::Fill(counts,0);last_cell=0;
    for(T_FACE face(local_grid,face_region);face;face++)if(!psi_N(face.Face())){
        CELL_ITERATOR& cell=cell_side?face.cell2:face.cell1;
        int c=cell.Cell(),color=filled_region_colors(cell.Cell());
        if(c!=last_cell && counts.Valid_Index(color)){last_cell=c;
            partitions(color).boundary_indices(side)(++counts(color))=cell_index_to_matrix_index(c);
            if(cell.Long()) partitions(color).boundary_indices(side)(++counts(color))=cell_index_to_matrix_index(c+1);}}
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID> void LAPLACE_RLE_MPI<T_GRID>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,ARRAY<int>&,ARRAY<ARRAY<int> >&,const VECTOR<T,2>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_RLE_MPI<T_GRID>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,ARRAY<int>&,ARRAY<ARRAY<int> >&,const VECTOR<T,3>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
#endif

//#####################################################################
template LAPLACE_RLE_MPI<RLE_GRID_2D<float> >::LAPLACE_RLE_MPI(LAPLACE_RLE<RLE_GRID_2D<float> >&);
template LAPLACE_RLE_MPI<RLE_GRID_3D<float> >::LAPLACE_RLE_MPI(LAPLACE_RLE<RLE_GRID_3D<float> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template LAPLACE_RLE_MPI<RLE_GRID_2D<double> >::LAPLACE_RLE_MPI(LAPLACE_RLE<RLE_GRID_2D<double> >&);
template LAPLACE_RLE_MPI<RLE_GRID_3D<double> >::LAPLACE_RLE_MPI(LAPLACE_RLE<RLE_GRID_3D<double> >&);
#endif
#endif
