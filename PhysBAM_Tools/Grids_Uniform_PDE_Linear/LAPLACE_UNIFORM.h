//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_UNIFORM
//#####################################################################
#ifndef __LAPLACE_UNIFORM__
#define __LAPLACE_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/LAPLACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_UNIFORM_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T> class PCG_SPARSE;
template<class T> class SPARSE_MATRIX_FLAT_NXN;
template<class T_GRID> class LAPLACE_UNIFORM_MPI;

template<class T_GRID>
class LAPLACE_UNIFORM:public LAPLACE<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL T_FLOOD_FILL;typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
public:
    typedef TV VECTOR_T;
    typedef T_GRID GRID_T;

    using LAPLACE<T>::tolerance;using LAPLACE<T>::number_of_regions;using LAPLACE<T>::solve_neumann_regions;

    T_GRID grid;
    T_ARRAYS_SCALAR& u;
    T_ARRAYS_SCALAR f; // f will be modified and reused as b in Ax=b for PCG
    PCG_SPARSE<T> pcg;
    PCG_SPARSE_THREADED<TV>* pcg_threaded;
    T_ARRAYS_INT filled_region_colors;
    ARRAY<bool> filled_region_touches_dirichlet;
    T_FACE_ARRAYS_BOOL psi_N;
    T_FACE_ARRAYS_SCALAR psi_R;
    T_ARRAYS_BOOL psi_D;
    TV_BOOL periodic_boundary;
    LAPLACE_UNIFORM_MPI<T_GRID>* laplace_mpi;
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    T_ARRAYS_BOOL* psi_D_save_for_sph;
    T_FACE_ARRAYS_BOOL* psi_N_save_for_sph;
    bool enforce_compatibility;
    bool solve_single_cell_neumann_regions;
    bool use_psi_R;
    THREAD_QUEUE* thread_queue;
#ifdef USE_PTHREADS
    pthread_barrier_t barr;
    pthread_mutex_t lock;
#endif
public:

    LAPLACE_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const bool initialize_grid,const bool enforce_compatibility_input,THREAD_QUEUE* thread_queue_input=0);
    virtual ~LAPLACE_UNIFORM();

    bool All_Cell_Faces_Neumann(const TV_INT& cell_index) const
    {for(int axis=1;axis<=T_GRID::dimension;axis++)if(!psi_N.Component(axis)(cell_index) || !psi_N.Component(axis)(cell_index+TV_INT::Axis_Vector(axis))) return false;
    return true;}

    bool Any_Cell_Faces_Neumann(const TV_INT& cell_index) const
    {for(int axis=1;axis<=T_GRID::dimension;axis++)if(psi_N.Component(axis)(cell_index) || psi_N.Component(axis)(cell_index+TV_INT::Axis_Vector(axis))) return true;
    return false;}

    bool Any_Neighbor_Dirichlet(const TV_INT& cell_index) const
    {for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
        if((!psi_N.Component(axis)(cell_index) && psi_D(cell_index-offset)) || (!psi_N.Component(axis)(cell_index+offset) && psi_D(cell_index+offset))) return true;}
    return false;}

//#####################################################################
    void Set_Neumann_Outer_Boundaries();
    void Set_Dirichlet_Outer_Boundaries();
    virtual void Initialize_Grid(const T_GRID& mac_grid_input);
    virtual void Solve(const T time=0,const bool solution_regions_already_computed=false);
    void Set_Threaded_Boundary(RANGE<TV_INT>& domain,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N);
    void Solve_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& threaded_grid,const T time=0);
    virtual void Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index);
    virtual void Find_A_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& row_counts);
    virtual void Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Find_Solution_Regions();
    void Compute_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Compute_Matrix_Indices(const RANGE<TV_INT>& domain,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Compute_Matrix_Indices_Threaded(ARRAY<RANGE<TV_INT> >& domains,ARRAY<ARRAY<INTERVAL<int> > >& interior_indices,ARRAY<ARRAY<ARRAY<INTERVAL<int> > > >& ghost_indices,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Solve_Subregion(ARRAY<TV_INT>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int color=0,ARRAY<int,TV_INT>* domain_index=0);
    void Solve_Subregion(ARRAY<INTERVAL<int> >& interior_indices,ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,ARRAY<TV_INT>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int color=0,ARRAY<int,TV_INT>* domain_index=0);
    void Build_Single_Solution_Region(T_ARRAYS_BOOL& solve);
    void Use_Psi_R();
//#####################################################################
};
}
#endif
