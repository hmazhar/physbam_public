//#####################################################################
// Copyright 2010 Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_THREADED
//#####################################################################
#ifndef __PCG_SPARSE_THREADED__
#define __PCG_SPARSE_THREADED__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/PTHREAD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class TV>
class PCG_SPARSE_THREADED:public PCG_SPARSE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef PCG_SPARSE<T> BASE;
public:
    THREAD_QUEUE& thread_queue;
    bool sum_int_needs_init,sum_needs_init,max_needs_init,barr_needs_init;
    bool needs_sync;
    T max_val;
    int sum_int,number_of_threads;
    VECTOR_ND<T> temp,p;
    ARRAY<T> sum_array;
#ifdef USE_PTHREADS
    pthread_mutex_t barr_lock,sum_lock,max_lock,sum_int_lock;
    pthread_mutex_t lock;
    pthread_barrier_t barr;
#endif

    using BASE::remove_null_space_solution_component;using BASE::enforce_compatibility;using BASE::maximum_iterations;
    using BASE::show_results;using BASE::show_residual;
    using BASE::incomplete_cholesky;using BASE::modified_incomplete_cholesky;using BASE::modified_incomplete_cholesky_coefficient;
    using BASE::preconditioner_zero_tolerance;using BASE::preconditioner_zero_replacement;

    PCG_SPARSE_THREADED(THREAD_QUEUE* thread_queue_input)
        :thread_queue(*thread_queue_input),sum_int_needs_init(false),sum_needs_init(false),max_needs_init(false),barr_needs_init(true),needs_sync(true),max_val(0),sum_int(0),sum_array(thread_queue.Number_Of_Threads())
    {
#ifdef USE_PTHREADS
        pthread_mutex_init(&sum_lock,0);
        pthread_mutex_init(&sum_int_lock,0);
        pthread_mutex_init(&max_lock,0);
        pthread_mutex_init(&barr_lock,0);
        pthread_mutex_init(&lock,0);
#endif
    }

    ~PCG_SPARSE_THREADED()
    {}

    void Init_Barriers()
    {
#ifdef USE_PTHREADS
        pthread_mutex_lock(&barr_lock);
        if(barr_needs_init){
            pthread_barrier_destroy(&barr);
            pthread_barrier_init(&barr,0,number_of_threads);
            barr_needs_init=false;}
        pthread_mutex_unlock(&barr_lock);
        pthread_barrier_wait(&barr);
        barr_needs_init=true;
#endif
    }

    int Global_Sum_Int(int input)
    {
#ifdef USE_PTHREADS
        pthread_mutex_lock(&sum_int_lock);
        if(sum_int_needs_init){sum_int=input;sum_int_needs_init=false;}
        else sum_int+=input;
        pthread_mutex_unlock(&sum_int_lock);
        pthread_barrier_wait(&barr);
        sum_int_needs_init=true;
        return sum_int;
#else
        return input;
#endif
    }


    T Global_Sum(T input,int tid)
    {
#ifdef USE_PTHREADS
        sum_array(tid)=input;
        pthread_barrier_wait(&barr);
        T sum=0;for(int i=1;i<=number_of_threads;i++) sum+=sum_array(i);
        return sum;
#else
        return input;
#endif
    }

    T Global_Max(T input)
    {
#ifdef USE_PTHREADS
        pthread_mutex_lock(&max_lock);
        if(max_needs_init){max_val=input;max_needs_init=false;}
        else max_val=max(max_val,input);
        pthread_mutex_unlock(&max_lock);
        pthread_barrier_wait(&barr);
        max_needs_init=true;
        return max_val;
#else
        return input;
#endif
    }

    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A_matrix,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance=1e-7,const bool recompute_preconditioner=true)
    {
        PHYSBAM_FATAL_ERROR();
    }

//#####################################################################
    void Solve(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance);
    void Solve_In_Parts(DOMAIN_ITERATOR_THREADED_ALPHA<PCG_SPARSE_THREADED<TV>,TV>& threaded_iterator,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance);
    void Solve_Part_One(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,VECTOR_ND<T>& x,VECTOR_ND<T>& b,ARRAY<VECTOR_ND<T> >& z_interior,ARRAY<VECTOR_ND<T> >& x_interior,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& p_interior,ARRAY<VECTOR_ND<T> >& temp_interior);
    void Solve_Part_Two(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& temp_interior);
    void Solve_Part_Three(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<SPARSE_MATRIX_FLAT_NXN<T>*>& C);
    void Solve_Part_Four(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& z_interior,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& temp_interior,ARRAY<SPARSE_MATRIX_FLAT_NXN<T>*>& C);
    void Solve_Part_Five(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& z_interior,ARRAY<VECTOR_ND<T> >& p_interior,T rho,T rho_old,int iteration);
    void Solve_Part_Six(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,const SPARSE_MATRIX_FLAT_NXN<T>& A);
    void Solve_Part_Seven(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& x_interior,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& p_interior,ARRAY<VECTOR_ND<T> >& temp_interior,T alpha);
    void Solve_Distribute(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior,const T sum);
    void Solve_Sum(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior,ARRAY<T>& sum);
    void Solve_Max(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior,ARRAY<T>& sum);
    void Solve_Dot(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior_1,ARRAY<VECTOR_ND<T> >& interior_2,ARRAY<T>& sum);
    void Solve_In_Parts(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance);
    void Threaded_Subtract(VECTOR_ND<T>& vector,const T sum,int start_index,int end_index);
    void Threaded_Sum(VECTOR_ND<T>& vector,ARRAY<T>& sum,int start_index,int end_index,int tid);
    void Threaded_Dot(VECTOR_ND<T>& vector1,VECTOR_ND<T>& vector2,ARRAY<T>& sum,int start_index,int end_index,int tid);
    void Threaded_Max(VECTOR_ND<T>& vector,ARRAY<T>& sum,int start_index,int end_index,int tid);
    void Threaded_Part_One(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T sum,int start_index,int end_index);
    void Threaded_Part_Two(VECTOR_ND<T>& z,T rho,T rho_old,int iteration,int start_index,int end_index);
    void Threaded_Part_Three(SPARSE_MATRIX_FLAT_NXN<T>& A,int start_index,int end_index);
    void Threaded_Part_Four(VECTOR_ND<T>& x,VECTOR_ND<T>& b,T alpha,int start_index,int end_index);
//#####################################################################
};
}
#endif
