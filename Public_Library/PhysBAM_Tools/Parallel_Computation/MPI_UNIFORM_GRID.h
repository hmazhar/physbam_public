//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_UNIFORM_GRID
//#####################################################################
#ifndef __MPI_UNIFORM_GRID__
#define __MPI_UNIFORM_GRID__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>
namespace PhysBAM{

template<class T_GRID>
class MPI_UNIFORM_GRID:public MPI_GRID<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<RANGE<TV_INT> >::TYPE T_ARRAYS_BOX_INT;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:
    typedef T_GRID GRID_T;

    typedef MPI_GRID<T_GRID> BASE;
    using BASE::Exchange_Boundary_Cell_Data;using BASE::local_grid;using BASE::global_grid;using BASE::comm;using BASE::all_neighbor_ranks;using BASE::side_neighbor_ranks;
    using BASE::Get_Unique_Tag;using BASE::Get_Send_Tag;using BASE::Get_Recv_Tag;using BASE::Wrap_Offset;using BASE::all_neighbor_directions;using BASE::all_coordinates;
    using BASE::Restrict_Grid;

    THREADED_UNIFORM_GRID<T_GRID>* threaded_grid;

    MPI_UNIFORM_GRID(T_GRID& local_grid_input,const int number_of_ghost_cells_input,const bool skip_initialization=false,const TV_INT& processes_per_dimension=TV_INT(),
        const TV_BOOL& periodic_input=TV_BOOL(),MPI::Group* group_input=0);
    
    MPI_UNIFORM_GRID(ARRAY<THREAD_PACKAGE>& buffers_input,const int tid_input,const int number_of_threads,T_GRID& local_grid_input,const int number_of_ghost_cells_input,
        const bool skip_mpi=true,const bool skip_initialization=false,const TV_INT& processes_per_dimension=TV_INT(),const TV_BOOL& periodic_input=TV_BOOL(),MPI::Group* group_input=0);

    template<class T_ARRAYS> void Exchange_Boundary_Cell_Data(T_ARRAYS& data,const int bandwidth,const bool include_corners=true) const
    {if(threaded_grid) threaded_grid->Exchange_Boundary_Cell_Data(data,bandwidth,include_corners);
    else MPI_GRID<T_GRID>::Exchange_Boundary_Cell_Data(*this,data,bandwidth,include_corners);}
    
    template<class T_ARRAYS> void Union_Common_Face_Data(T_ARRAYS& data) const
    {MPI_GRID<T_GRID>::Union_Common_Face_Data(*this,data);}

    template<class T_FACE_ARRAYS_2>
    void Exchange_Boundary_Face_Data(T_FACE_ARRAYS_2& data,const int bandwidth) const
    {if(threaded_grid) threaded_grid->Exchange_Boundary_Face_Data(data,bandwidth);
    else MPI_GRID<T_GRID>::Exchange_Boundary_Face_Data(*this,data,bandwidth);}

    template<class T_FACE_ARRAYS_2> void Average_Common_Face_Data(T_FACE_ARRAYS_2& data) const
    {if(threaded_grid) threaded_grid->Average_Common_Face_Data(data);
    else MPI_GRID<T_GRID>::Average_Common_Face_Data(*this,data);}

    template<class T_FACE_ARRAYS_2> void Copy_Common_Face_Data(T_FACE_ARRAYS_2& data) const
    {MPI_GRID<T_GRID>::Copy_Common_Face_Data(*this,data);}

    template<class T_FACE_ARRAYS_2> void Assert_Common_Face_Data(T_FACE_ARRAYS_2& data,const T tolerance=0) const
    {if(threaded_grid) threaded_grid->Assert_Common_Face_Data(data,tolerance);
    else MPI_GRID<T_GRID>::Assert_Common_Face_Data(*this,data,tolerance);}

    RANGE<TV_INT> Face_Sentinels(const int axis) const
    {return RANGE<TV_INT>(TV_INT(),TV_INT::Axis_Vector(axis));}

    RANGE<TV_INT> Parallel_Face_Sentinels(const int axis) const
    {return Face_Sentinels(axis);}
    
    void Synchronize_Dt(T& dt) const
    {if(threaded_grid) threaded_grid->Synchronize_Dt(dt);
    else MPI_GRID<T_GRID>::Synchronize_Dt(dt);}
    
    void Initialize(VECTOR<VECTOR<bool,2>,T_GRID::dimension>& domain_walls)
    {if(threaded_grid) threaded_grid->Initialize(domain_walls);
    else BASE::Initialize(domain_walls);}
    
    bool Neighbor(const int axis,const int axis_side) const
    {if(threaded_grid) return threaded_grid->Neighbor(axis,axis_side);
    else return BASE::Neighbor(axis,axis_side);}

//#####################################################################
    T_GRID Get_Non_Overlapping_Face_Grid(const int axis) const;
    template<class T_ARRAYS> bool Gather_Cell_Data(const T_ARRAYS& local_data,T_ARRAYS& global_data) const;
    template<class T_ARRAYS> void Scatter_Cell_Data(const T_ARRAYS& global_data,T_ARRAYS& local_data) const;
    template<class T2> MPI_PACKAGE Package_Cell_Data(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const RANGE<TV_INT>& region) const;
    template<class T_FACE_ARRAYS2> MPI_PACKAGE Package_Face_Data(T_FACE_ARRAYS2& data,const ARRAY<RANGE<TV_INT> >& region) const;
    template<class T_FACE_ARRAYS2> MPI_PACKAGE Package_Common_Face_Data(T_FACE_ARRAYS2& data,const int axis,const RANGE<TV_INT>& region) const;
//#####################################################################
};
}
#endif
