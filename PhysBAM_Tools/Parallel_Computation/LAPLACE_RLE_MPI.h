//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_RLE_MPI
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __LAPLACE_RLE_MPI__
#define __LAPLACE_RLE_MPI__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_MPI.h>
namespace PhysBAM{

class GRAPH;
template<class T_GRID> class LAPLACE_COLLIDABLE_RLE;

template<class T_GRID>
class LAPLACE_RLE_MPI:public LAPLACE_MPI<T_GRID>
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef MPI_RLE_GRID<T_GRID> T_MPI_GRID;
    typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename T_GRID::HORIZONTAL_GRID T_HORIZONTAL_GRID;

    typedef LAPLACE_MPI<T_GRID> BASE;
    using BASE::mpi_grid;using BASE::psi_N;using BASE::filled_region_colors;using BASE::local_grid;using BASE::number_of_regions;using BASE::filled_region_ranks;using BASE::partitions;

    LAPLACE_RLE_MPI(LAPLACE_RLE<T_GRID>& laplace)
        :LAPLACE_MPI<T_GRID>(laplace)
    {}

    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& matrix_index_to_cell_index_array) PHYSBAM_OVERRIDE
    {Find_Matrix_Indices(filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,TV());}

//#####################################################################
private:
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&);
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&);
    void Find_Matrix_Indices_In_Region(const int region_index,const typename T_GRID::BOX_HORIZONTAL_INT& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,
        ARRAY<ARRAY<int> >& matrix_index_to_cell_index_array);
    template<class T_FACE> void Find_Boundary_Indices_In_Region(const int side,const typename T_GRID::BOX_HORIZONTAL_INT& region,ARRAY<int>& cell_index_to_matrix_index);
//#####################################################################
};
}
#endif
#endif
