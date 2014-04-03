//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_GRID
//#####################################################################
#ifndef __MPI_GRID__
#define __MPI_GRID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;
template<class TV> class RANGE;
template<class TV> class GRID;

template<class T_GRID>
class MPI_GRID:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
public:
    T_GRID& local_grid;
    int number_of_ghost_cells;
    T_GRID global_grid;
    int rank;
    int number_of_processes;
    MPI::Intracomm* comm;
    MPI::Group* group;
    T_GRID process_grid;
    T_ARRAYS_INT process_ranks;
    TV_INT coordinates;
    ARRAY<TV_INT> all_coordinates;
    ARRAY<int> all_neighbor_ranks,side_neighbor_ranks; // all_neighbor_ranks includes all 1 ring neighbors
    ARRAY<TV_INT> all_neighbor_directions,side_neighbor_directions;
    TV_INT local_to_global_offset;
    T_ARRAYS_INT local_cell_index_to_global_column_index_map;
    ARRAY<ARRAY<int> > boundaries;
private:
    mutable int current_tag;
public:
    TV_BOOL periodic;
    bool ignore_boundary_faces;

    MPI_GRID(T_GRID& local_grid_input,const int number_of_ghost_cells_input,const bool skip_initialization=false,const TV_INT& processes_per_dimension=TV_INT(),
        const TV_BOOL& periodic_input=TV_BOOL(),MPI::Group* group_input=0);
    ~MPI_GRID();

    int Number_Of_Processors() const
    {return number_of_processes;}

    int Get_Unique_Tag() const
    {return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

    int Get_Send_Tag(const TV_INT& direction) const
    {STATIC_ASSERT(TV_INT::m<=3);int tag=0;
    for(int i=1;i<=direction.m;i++){assert(abs(direction[i])<=1);tag=3*tag+direction[i]+1;}
    return tag;}

    int Get_Recv_Tag(const TV_INT& direction) const
    {return Get_Send_Tag(-direction);}

protected:
    TV Wrap_Offset(const TV_INT& direction) const // offset to add to translate into space of adjacent processor
    {TV offset;TV_INT neighbor_coordinates=coordinates+direction;
    for(int axis=1;axis<=offset.m;axis++)if(periodic[axis] && (neighbor_coordinates[axis]<1 || neighbor_coordinates[axis]>process_grid.counts[axis]))
        offset[axis]=-direction[axis]*global_grid.domain.Edge_Lengths()[axis];
    return offset;}
public:

//#####################################################################
    void Initialize(VECTOR<VECTOR<bool,2>,T_GRID::dimension>& domain_walls);
    bool Neighbor(const int axis,const int axis_side) const;
    void Split_Grid(const TV_INT& processes_per_dimension);
    T_GRID Restrict_Grid(const TV_INT& coordinates) const;
    void Synchronize_Dt(T& dt) const;
    void Synchronize_J_Bounds(int& jmin,int& jmax) const;
    void Sync_Common_Face_Weights_To(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >,FACE_INDEX<TV::dimension> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >,FACE_INDEX<TV::dimension> >& weights_from,const int ghost_cells);
    void Sync_Common_Face_Weights_From(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >,FACE_INDEX<TV::dimension> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >,FACE_INDEX<TV::dimension> >& weights_from,const int ghost_cells);
    void Sync_Common_Cell_Weights_To(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells);
    void Sync_Common_Cell_Weights_From(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells);
    template<class T_MPI_GRID,class T2> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const int bandwidth,const bool include_corners=true) const;
    template<class T_MPI_GRID,class T2> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const T_GRID& local_grid,ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const int bandwidth,
        const bool include_corners=true) const;
    template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,const int bandwidth,const bool include_corners=true) const;
    template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const T_GRID& local_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,const int bandwidth,
        const bool include_corners=true) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Exchange_Boundary_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data,const int bandwidth) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Average_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Copy_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Assert_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data,const T tolerance=0) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS_BOOL> void Union_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS_BOOL& data) const;
    void Find_Boundary_Regions(ARRAY<RANGE<TV_INT> >& regions,const RANGE<TV_INT>& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions=true) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions,const GRID<VECTOR<T,1> >& local_grid) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions,const GRID<VECTOR<T,2> >& local_grid) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,3> > >& regions,const RANGE<VECTOR<int,3> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions,const GRID<VECTOR<T,3> >& local_grid) const;
    RANGE<TV_INT> Find_Region_Box(const int processor,const RANGE<TV_INT>& sentinels,const int band) const;
    template<class T2> void Reduce_Add(const T2& input,T2& output) const;
    template<class T2> T2 Reduce_Add(const T2& local_value) const;

protected:
    VECTOR<int,1> Split_Grid(const GRID<VECTOR<T,1> >& global_grid,const VECTOR<int,1>& processes_per_dimension);
    VECTOR<int,2> Split_Grid(const GRID<VECTOR<T,2> >& global_grid,const VECTOR<int,2>& processes_per_dimension);
    VECTOR<int,3> Split_Grid(const GRID<VECTOR<T,3> >& global_grid,const VECTOR<int,3>& processes_per_dimension);
    static void Split_Dimension(const int m,const int processes,ARRAY<int>& boundaries);
    void Initialize_Communicator(const bool manual,MPI::Group* group);
//#####################################################################
};
}
#endif
