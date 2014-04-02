//#####################################################################
// Copyright 2010 Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_MPI_THREADED
//#####################################################################
#ifndef __PCG_SPARSE_MPI_THREADED__
#define __PCG_SPARSE_MPI_THREADED__
#ifdef USE_MPI
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
namespace PhysBAM{

template<class TV>
class PCG_SPARSE_MPI_THREADED:public PCG_SPARSE_MPI<GRID<VECTOR<typename TV::SCALAR,TV::dimension> > >
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef PCG_SPARSE_MPI<GRID<TV> > BASE;
public:
    using BASE::pcg;using BASE::partition;using BASE::Fill_Ghost_Cells;

    PCG_SPARSE_THREADED<TV>& pcg_threaded;

    PCG_SPARSE_MPI_THREADED(PCG_SPARSE_THREADED<TV>& pcg_threaded_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input)
        :BASE(pcg_threaded_input,comm_input,partition_input),pcg_threaded(pcg_threaded_input)
    {}

    ~PCG_SPARSE_MPI_THREADED()
    {}

    int Global_Sum_Int(int input)
    {int local_sum=pcg_threaded.Global_Sum_Int(input);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&pcg_threaded.sum_int_lock);
    if(pcg_threaded.needs_sync){pcg_threaded.sum_int=BASE::Global_Sum(local_sum);pcg_threaded.needs_sync=false;}
    pthread_mutex_unlock(&pcg_threaded.sum_int_lock);
    pthread_barrier_wait(&pcg_threaded.barr);
    pcg_threaded.needs_sync=true;
    pthread_barrier_wait(&pcg_threaded.barr);
    return pcg_threaded.sum_int;
#else
    return BASE::Global_Sum(local_sum);
#endif
    }

    T Global_Sum(T input,int tid)
    {T local_sum=pcg_threaded.Global_Sum(input,tid);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&pcg_threaded.sum_lock);
    T sum=0;if(pcg_threaded.needs_sync){sum=BASE::Global_Sum(local_sum);pcg_threaded.needs_sync=false;}
    pthread_mutex_unlock(&pcg_threaded.sum_lock);
    pthread_barrier_wait(&pcg_threaded.barr);
    pcg_threaded.needs_sync=true;
    pthread_barrier_wait(&pcg_threaded.barr);
    return sum;
#else
    return BASE::Global_Sum(local_sum);
#endif
    }

    T Global_Max(T input)
    {T local_max=pcg_threaded.Global_Max(input);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&pcg_threaded.max_lock);
    if(pcg_threaded.needs_sync){pcg_threaded.max_val=BASE::Global_Max(local_max);pcg_threaded.needs_sync=false;}
    pthread_mutex_unlock(&pcg_threaded.max_lock);
    pthread_barrier_wait(&pcg_threaded.barr);
    pcg_threaded.needs_sync=true;
    pthread_barrier_wait(&pcg_threaded.barr);
    return pcg_threaded.max_val;
#else
    return BASE::Global_Max(local_max);
#endif
    }

    void Fill_Ghost_Cells(VECTOR_ND<T>& v)
    {
#ifdef USE_PTHREADS
    pthread_barrier_wait(&pcg_threaded.barr);
    pthread_mutex_lock(&pcg_threaded.lock);
    if(pcg_threaded.needs_sync){BASE::Fill_Ghost_Cells(v);pcg_threaded.needs_sync=false;}
    pthread_mutex_unlock(&pcg_threaded.lock);
    pthread_barrier_wait(&pcg_threaded.barr);
    pcg_threaded.needs_sync=true;
    pthread_barrier_wait(&pcg_threaded.barr);
#else
    BASE::Fill_Ghost_Cells(v);
#endif
    }

    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A_matrix,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance=1e-7,const bool recompute_preconditioner=true)
    {
        PHYSBAM_FATAL_ERROR();
    }

//#####################################################################
    void Solve(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance);
//#####################################################################
};
}
#endif
#endif
