//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_MPI
//#####################################################################
#ifndef __PCG_SPARSE_MPI__
#define __PCG_SPARSE_MPI__

#ifdef USE_MPI

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>
namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;

template<class T_GRID>
class PCG_SPARSE_MPI:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
public:
    PCG_SPARSE<T>& pcg;
    MPI::Intracomm& comm;
    THREADED_UNIFORM_GRID<T_GRID>* thread_grid;
    SPARSE_MATRIX_PARTITION& partition;
    ARRAY<MPI::Datatype> boundary_datatypes,ghost_datatypes;
    ARRAY<ARRAY<int> > columns_to_send;
    ARRAY<ARRAY<int> > columns_to_receive;

    PCG_SPARSE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input)
        :pcg(pcg_input),comm(comm_input),partition(partition_input)
    {}

    ~PCG_SPARSE_MPI()
    {MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);}

    template<class TYPE> TYPE Global_Sum(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::SUM,comm);return output;}

    template<class TYPE> TYPE Global_Max(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::MAX,comm);return output;}

    void Fill_Ghost_Cells(VECTOR_ND<T>& v)
    {ARRAY<MPI::Request> requests;requests.Preallocate(2*partition.number_of_sides);
    for(int s=1;s<=partition.number_of_sides;s++)if(boundary_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append(comm.Isend(v.x-1,1,boundary_datatypes(s),partition.neighbor_ranks(s),s));
    for(int s=1;s<=partition.number_of_sides;s++)if(ghost_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append(comm.Irecv(v.x-1,1,ghost_datatypes(s),partition.neighbor_ranks(s),((s-1)^1)+1));
    MPI_UTILITIES::Wait_All(requests);}
    
//#####################################################################
    void Serial_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const int tag,const T tolerance=1e-7);
    void Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x_local,VECTOR_ND<T>& b_local,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries,
        const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Find_Ghost_Regions(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries);
    void Find_Ghost_Regions_Threaded(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries);
    void Fill_Ghost_Cells_Far(VECTOR_ND<T>& x);
    void Fill_Ghost_Cells_Threaded(VECTOR_ND<T>& x);    
    void Initialize_Datatypes();
//#####################################################################
};
}
#endif
#endif
