//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THREADED_UNIFORM_GRID
//#####################################################################
#ifndef __THREADED_UNIFORM_GRID__
#define __THREADED_UNIFORM_GRID__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/PTHREAD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_PACKAGE.h>

//TODO: Merge this with MPI_UNIFORM_GRID
namespace PhysBAM{

template<class T_GRID>
class THREADED_UNIFORM_GRID:public MPI_GRID<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;
    typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<RANGE<TV_INT> >::TYPE T_ARRAYS_BOX_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
public:
    typedef T_GRID GRID_T;

    typedef MPI_GRID<T_GRID> BASE;
    using BASE::number_of_processes;using BASE::global_grid;using BASE::local_grid;using BASE::coordinates;using BASE::boundaries;using BASE::periodic;
    using BASE::side_neighbor_ranks;using BASE::side_neighbor_directions;using BASE::all_neighbor_ranks;using BASE::all_neighbor_directions;using BASE::rank;
    using BASE::local_to_global_offset;using BASE::all_coordinates;using BASE::local_cell_index_to_global_column_index_map;using BASE::process_ranks;using BASE::process_grid;
    using BASE::Restrict_Grid;using BASE::Find_Boundary_Regions;

    int tid;
    ARRAY<VECTOR<int,2> > global_column_index_boundaries;    
    ARRAY<THREAD_PACKAGE>& buffers; //global memory
#ifdef USE_PTHREADS
    mutable pthread_mutex_t* lock;
    mutable pthread_barrier_t* barr;
#endif

    THREADED_UNIFORM_GRID(ARRAY<THREAD_PACKAGE>& buffers_input,const int tid_input,const int number_of_threads,T_GRID& local_grid_input,const int number_of_ghost_cells_input,
        const bool skip_initialization=false,const TV_INT& processes_per_dimension=TV_INT(),const TV_BOOL& periodic_input=TV_BOOL());

    RANGE<TV_INT> Face_Sentinels(const int axis) const
    {return RANGE<TV_INT>(TV_INT(),TV_INT::Axis_Vector(axis));}

    bool Neighbor(const int axis,const int axis_side) const
    {return side_neighbor_ranks(2*(axis-1)+axis_side)!=-1;}
    
    template<class T2> THREAD_PACKAGE Package_Cell_Data(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const RANGE<TV_INT>& send_region) const
    {
        int size=0;for(CELL_ITERATOR iterator(local_grid,send_region);iterator.Valid();iterator.Next()) size+=data.Pack_Size();
        int position=0;THREAD_PACKAGE pack(size);pack.send_tid=rank;
        for(CELL_ITERATOR iterator(local_grid,send_region);iterator.Valid();iterator.Next()) data.Pack(pack.buffer,position,iterator.Cell_Index());
        return pack;
    }

    template<class T2> void Sync_Particles(const ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& local_data,ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& global_data) const
    {
        RANGE<TV_INT> domain=local_grid.Domain_Indices();
        for(int axis=1;axis<=TV::dimension;axis++) if(domain.max_corner(axis)+local_to_global_offset(axis)==global_grid.Domain_Indices().max_corner(axis)) domain.max_corner(axis)++;
        for(NODE_ITERATOR iterator(local_grid,domain);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();global_data(node+local_to_global_offset)=local_data(node);}
    }

    void Initialize(VECTOR<VECTOR<bool,2>,T_GRID::dimension>& domain_walls);
    void Synchronize_Dt(T& dt) const;
    void All_Reduce(bool& flag) const;
    void Allgather(ARRAY<int>& data) const;
    template<class T2> void Exchange_Boundary_Cell_Data(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const int bandwidth,const bool include_corners=true) const;
    template<class T2> void Exchange_Boundary_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data,const int bandwidth) const;
    template<class T2> void Average_Common_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data) const;
    template<class T2> void Assert_Common_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data,const T tolerance) const;
    template<class T2> void Sync_Scalar(const ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& local_data,ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& global_data) const;
    template<class T2> void Sync_Face_Scalar(const ARRAY<T2,FACE_INDEX<TV::dimension> >& local_data,ARRAY<T2,FACE_INDEX<TV::dimension> >& global_data) const;
    template<class T2> void Distribute_Scalar(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& local_data,const ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& global_data) const;
    template<class T2> void Distribute_Face_Scalar(ARRAY<T2,FACE_INDEX<TV::dimension> >& local_data,const ARRAY<T2,FACE_INDEX<TV::dimension> >& global_data) const;
};
}
#endif
