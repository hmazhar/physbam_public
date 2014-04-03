//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class POLICY> RLE_GRID<T,POLICY>::
RLE_GRID()
    :negative_bandwidth(0),positive_bandwidth(0),number_of_cells(0),number_of_faces(0),number_of_blocks(0),mpi_grid(0),jmin(0),jmax(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class POLICY> RLE_GRID<T,POLICY>::
~RLE_GRID()
{}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T,class POLICY> void RLE_GRID<T,POLICY>::
Clean_Memory()
{
    columns.Clean_Memory();Topology_Changed();
}
//#####################################################################
// Function Transfer
//#####################################################################
template<class T,class POLICY> void RLE_GRID<T,POLICY>::
Transfer(RLE_GRID<T,POLICY>& old_grid,RLE_GRID<T,POLICY>& new_grid)
{
    Transfer_Noncolumn_Data(old_grid,new_grid);
    T_ARRAYS_HORIZONTAL_ARRAY_RUN::Exchange_Arrays(old_grid.columns,new_grid.columns);
    T_ARRAYS_HORIZONTAL_INT::Exchange_Arrays(old_grid.block_run_offsets,new_grid.block_run_offsets);
    ARRAY<T_BLOCK_RUN>::Exchange_Arrays(old_grid.block_runs,new_grid.block_runs);
    old_grid.Clean_Memory();
} 
//#####################################################################
// Function Transfer_Noncolumn_Data
//#####################################################################
template<class T,class POLICY> void RLE_GRID<T,POLICY>::
Transfer_Noncolumn_Data(const RLE_GRID<T,POLICY>& old_grid,RLE_GRID<T,POLICY>& new_grid)
{
    new_grid.Clean_Memory();
    new_grid.uniform_grid=old_grid.uniform_grid;new_grid.horizontal_grid=old_grid.horizontal_grid;new_grid.number_of_ghost_cells=old_grid.number_of_ghost_cells;
    new_grid.minimum_vertical_space=old_grid.minimum_vertical_space;new_grid.minimum_long_run_length=old_grid.minimum_long_run_length;
    new_grid.long_run_cells=old_grid.long_run_cells;new_grid.long_run_faces_horizontal=old_grid.long_run_faces_horizontal;
    new_grid.negative_bandwidth=old_grid.negative_bandwidth;new_grid.positive_bandwidth=old_grid.positive_bandwidth;
    new_grid.number_of_cells=old_grid.number_of_cells;new_grid.number_of_faces=old_grid.number_of_faces;new_grid.last_face_numbers=old_grid.last_face_numbers;
    new_grid.number_of_blocks=old_grid.number_of_blocks;new_grid.jmin=old_grid.jmin;new_grid.jmax=old_grid.jmax;new_grid.domain=old_grid.domain;
    new_grid.mpi_grid=old_grid.mpi_grid;
} 
//#####################################################################
// Function Topology_Changed
//#####################################################################
template<class T,class POLICY> void RLE_GRID<T,POLICY>::
Topology_Changed()
{
    block_run_offsets.Clean_Memory();block_runs.Clean_Memory();short_cell_neighbors.Clean_Memory();short_face_neighbors.Clean_Memory();
    short_face_locations.Clean_Memory();
}
//#####################################################################
// Function Union
//#####################################################################
template<class T,class POLICY> void RLE_GRID<T,POLICY>::
Union(const T_GRID& grid1,const T_GRID& grid2,T_GRID& union_grid,const T_ARRAYS_HORIZONTAL_INT& ground_j)
{
    if(grid1.domain!=grid2.domain) PHYSBAM_FATAL_ERROR();
    Transfer_Noncolumn_Data(grid1,union_grid);
    union_grid.columns.Resize(union_grid.horizontal_grid.Get_MAC_Grid().Domain_Indices(union_grid.number_of_ghost_cells+1));
    for(HORIZONTAL_CELL_ITERATOR iterator(union_grid.horizontal_grid,union_grid.number_of_ghost_cells+1);iterator.Valid();iterator.Next()){
        ARRAY<RLE_RUN> column_union;
        TV_HORIZONTAL_INT column=iterator.Cell_Index();
        RLE_RUN::Column_Union(grid1.columns(column),grid2.columns(column),column_union,union_grid.minimum_long_run_length);
        RLE_RUN::Split_At_Ground(column_union,union_grid.columns(column),ground_j(column),union_grid.minimum_long_run_length);}
    union_grid.Compute_Auxiliary_Information();
}
//#####################################################################
// Function Put_Ghost_Helper
//#####################################################################
namespace{
template<class T2> struct Put_Ghost_Helper{template<class T_ITERATOR,class T_GRID> static void Apply(const T_GRID& grid,const T2& constant,ARRAY<T2>& array)
{
    ARRAY<RANGE<VECTOR<int,T_GRID::VECTOR_T::m-1> > > ghost_regions;grid.Find_Ghost_Regions(ghost_regions,T_ITERATOR::Sentinels(),grid.number_of_ghost_cells);
    for(int region=1;region<=ghost_regions.m;region++)for(RLE_GRID_SIMPLE_ITERATOR<T_GRID,T_ITERATOR> iterator(grid,ghost_regions(region));iterator;iterator++)array(iterator.index)=constant;
}};}
//#####################################################################
// Function Put_Ghost_Cells
//#####################################################################
template<class T,class POLICY> template<class T2> void RLE_GRID<T,POLICY>::
Put_Ghost_Cells(const T2& constant,ARRAY<T2>& array) const
{
    Put_Ghost_Helper<T2>::template Apply<CELL_ITERATOR>(static_cast<const T_GRID&>(*this),constant,array);
}
//#####################################################################
// Function Put_Ghost_Faces
//#####################################################################
template<class T,class POLICY> template<class T2> void RLE_GRID<T,POLICY>::
Put_Ghost_Faces(const T2& constant,ARRAY<T2>& array) const
{
    this->template Face_Loop<Put_Ghost_Helper<T2> >(static_cast<const T_GRID&>(*this),constant,array);
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class RLE_GRID<T,RLE_POLICY_##d##D<T> >; \
    template void RLE_GRID<T,RLE_POLICY_##d##D<T> >::Put_Ghost_Cells(const T&,ARRAY<T>&) const; \
    template void RLE_GRID<T,RLE_POLICY_##d##D<T> >::Put_Ghost_Cells(const bool&,ARRAY<bool>&) const; \
    template void RLE_GRID<T,RLE_POLICY_##d##D<T> >::Put_Ghost_Faces(const T&,ARRAY<T>&) const; \
    template void RLE_GRID<T,RLE_POLICY_##d##D<T> >::Put_Ghost_Faces(const bool&,ARRAY<bool>&) const;
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif
#endif
