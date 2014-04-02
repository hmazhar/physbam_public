//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_UNIFORM_MPI
//#####################################################################
#ifndef __LAPLACE_UNIFORM_MPI__
#define __LAPLACE_UNIFORM_MPI__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/PTHREAD.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
namespace PhysBAM{

class GRAPH;
template<class TV> class GRID;
template<class T_GRID> class LAPLACE_UNIFORM;

template<class T_GRID>
class LAPLACE_UNIFORM_MPI:public LAPLACE_MPI<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
public:
    typedef LAPLACE_MPI<T_GRID> BASE;
    using BASE::mpi_grid;using BASE::local_grid;using BASE::filled_region_ranks;using BASE::partitions;using BASE::number_of_regions;using BASE::solve_neumann_regions;using BASE::psi_N;
    using BASE::filled_region_touches_dirichlet;

    bool sum_int_needs_init;
    int sum_int;
    pthread_mutex_t sum_int_lock,lock;
    pthread_barrier_t barr;


    LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<T_GRID>& laplace)
        :LAPLACE_MPI<T_GRID>(laplace),sum_int_needs_init(false),sum_int(0)
    {pthread_mutex_init(&sum_int_lock,0);pthread_mutex_init(&lock,0);}

    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array) PHYSBAM_OVERRIDE
    {Find_Matrix_Indices(filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,TV());}

    int Global_Sum(int input)
    {
#ifdef USE_PTHREADS
        pthread_mutex_lock(&sum_int_lock);
        if(sum_int_needs_init){sum_int=input;sum_int_needs_init=false;}
        else sum_int+=input;
        pthread_mutex_unlock(&sum_int_lock);
        pthread_barrier_wait(&barr);
        sum_int_needs_init=true;
        pthread_barrier_wait(&barr);
        return sum_int;
#else
        return input;
#endif
    }

//#####################################################################
    void Find_Matrix_Indices_Threaded(ARRAY<RANGE<TV_INT> >& domains,ARRAY<ARRAY<INTERVAL<int> > >& interior_indices,ARRAY<ARRAY<ARRAY<INTERVAL<int> > > >& ghost_indices,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,LAPLACE_UNIFORM<T_GRID>* laplace);
private:
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,1>&);
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&);
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&);
    void Find_Matrix_Indices_In_Region(const int region_index,const RANGE<TV_INT>& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,
        ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array);
    void Find_Boundary_Indices_In_Region(const int side,const RANGE<TV_INT>& region,T_ARRAYS_INT& cell_index_to_matrix_index);
//#####################################################################
};
}
#endif
