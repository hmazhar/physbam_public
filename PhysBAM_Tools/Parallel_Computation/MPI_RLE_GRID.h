//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_RLE_GRID
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __MPI_RLE_GRID__
#define __MPI_RLE_GRID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>
namespace PhysBAM{

class RLE_RUN;
class MPI_PACKAGE;
template<class TV> class RANGE;
template<class T1,class T2> class PAIR;

template<class T_GRID>
class MPI_RLE_GRID:public MPI_GRID<typename T_GRID::HORIZONTAL_GRID>
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::RUN T_RUN;
    typedef typename T_GRID::ARRAYS_HORIZONTAL::template REBIND<ARRAY<RLE_RUN> >::TYPE T_ARRAYS_HORIZONTAL_ARRAY_RUN;
    typedef typename T_GRID::HORIZONTAL_GRID T_HORIZONTAL_GRID;
    typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;typedef typename T_GRID::BOX_HORIZONTAL T_BOX_HORIZONTAL;
    typedef typename GRID<TV>::VECTOR_INT TV_INT;typedef typename T_GRID::BLOCK_ITERATOR BLOCK_ITERATOR;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;

    typedef T_GRID GRID_T;

    typedef MPI_GRID<T_HORIZONTAL_GRID> BASE;
    using BASE::comm;using BASE::all_neighbor_ranks;using BASE::all_neighbor_directions;using BASE::Get_Unique_Tag;using BASE::Get_Send_Tag;using BASE::Get_Recv_Tag;
    using BASE::Find_Boundary_Regions;using BASE::rank;using BASE::Wrap_Offset;

    T_GRID& local_grid;
    GRID<TV> global_uniform_grid;

    MPI_RLE_GRID(T_GRID& local_grid_input,const bool skip_initialization=false,const TV_HORIZONTAL_INT& processes_per_dimension=TV_HORIZONTAL_INT(),const TV_BOOL& periodic_input=TV_BOOL())
        :MPI_GRID<T_HORIZONTAL_GRID>(local_grid_input.horizontal_grid,local_grid_input.number_of_ghost_cells,skip_initialization,processes_per_dimension,periodic_input.Horizontal_Vector()),
        local_grid(local_grid_input)
    {
        PHYSBAM_ASSERT(!periodic_input.y);
        if(skip_initialization) return;
        global_uniform_grid=local_grid.uniform_grid;
        Update_Local_Grid(MPI_GRID<T_HORIZONTAL_GRID>::local_grid);
    }

    virtual ~MPI_RLE_GRID()
    {}

private:
    static GRID<VECTOR<T,2> > Vertical_Grid(const GRID<VECTOR<T,2> >& vertical,const GRID<VECTOR<T,1> >& horizontal)
    {return GRID<VECTOR<T,2> >(horizontal.counts.x,vertical.counts.y,RANGE<VECTOR<T,2> >(horizontal.domain.min_corner.x,horizontal.domain.max_corner.x,vertical.domain.min_corner.y,vertical.domain.max_corner.y));}

    static GRID<VECTOR<T,3> > Vertical_Grid(const GRID<VECTOR<T,3> >& vertical,const GRID<VECTOR<T,2> >& horizontal)
    {return GRID<VECTOR<T,3> >(horizontal.counts.x,vertical.counts.y,horizontal.counts.y,RANGE<VECTOR<T,3> >(horizontal.domain.min_corner.x,horizontal.domain.max_corner.x,vertical.domain.min_corner.y,vertical.domain.max_corner.y,horizontal.domain.min_corner.y,horizontal.domain.max_corner.y));}
public:

    void Update_Local_Grid(const T_HORIZONTAL_GRID& local_horizontal_grid)
    {local_grid.Set_Uniform_Grid(Vertical_Grid(global_uniform_grid,local_horizontal_grid.Get_Regular_Grid()));}

    void Initialize(VECTOR<VECTOR<bool,2>,T_GRID::dimension-1>& horizontal_domain_walls)
    {MPI_GRID<T_HORIZONTAL_GRID>::Initialize(horizontal_domain_walls);}

    template<class T2> void Exchange_Boundary_Cell_Data(ARRAY_BASE<T2,ARRAY<T2> >& data,const int bandwidth,const bool include_corners=true) const
    {MPI_GRID<T_HORIZONTAL_GRID>::Exchange_Boundary_Cell_Data(*this,data,bandwidth,include_corners);}

    void Exchange_Boundary_Face_Data(ARRAY<T>& data,const int bandwidth) const
    {MPI_GRID<T_HORIZONTAL_GRID>::Exchange_Boundary_Face_Data(*this,data,bandwidth);}

    void Average_Common_Face_Data(ARRAY<T>& data) const
    {MPI_GRID<T_HORIZONTAL_GRID>::Average_Common_Face_Data(*this,data);}

    T_BOX_HORIZONTAL_INT Parallel_Face_Sentinels(const int horizontal_axis) const
    {return T_BOX_HORIZONTAL_INT(TV_HORIZONTAL_INT(),TV_HORIZONTAL_INT::Axis_Vector(horizontal_axis));}

    T_BOX_HORIZONTAL_INT Face_Sentinels(const int axis) const
    {return axis==2?T_BOX_HORIZONTAL_INT::Zero_Box():Parallel_Face_Sentinels(axis/2+1);}

    void Find_Boundary_Regions(ARRAY<T_BOX_HORIZONTAL_INT>& regions,const T_BOX_HORIZONTAL_INT& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions,const T_GRID& local_grid) const
    {Find_Boundary_Regions(regions,sentinels,skip_common_boundary,band,include_corners,include_ghost_regions);} // TODO: clean up the need for this forwarding function

//#####################################################################
    template<class T_ARRAYS_HORIZONTAL_COLUMN> void Exchange_Boundary_Columns(const int bandwidth,T_ARRAYS_HORIZONTAL_COLUMN& columns);
    template<class TS> MPI_PACKAGE Package_Cell_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const RANGE<VECTOR<int,1> >& region) const;
    template<class TS> MPI_PACKAGE Package_Cell_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const RANGE<VECTOR<int,2> >& region) const;
    template<class TS> MPI_PACKAGE Package_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const ARRAY<RANGE<VECTOR<int,1> > >& regions) const;
    template<class TS> MPI_PACKAGE Package_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const ARRAY<RANGE<VECTOR<int,2> > >& regions) const;
    template<class TS> MPI_PACKAGE Package_Common_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const int horizontal_axis,const RANGE<VECTOR<int,1> >& region) const;
    template<class TS> MPI_PACKAGE Package_Common_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const int horizontal_axis,const RANGE<VECTOR<int,2> >& region) const;
protected:
    template<class T_ARRAYS_HORIZONTAL_COLUMN> MPI::Request ISend_Columns(const T_ARRAYS_HORIZONTAL_COLUMN& columns,const ARRAY<T_BOX_HORIZONTAL_INT>& regions,const int neighbor,const int tag,
        ARRAY<char>& buffer) const;
    template<class T_ARRAYS_HORIZONTAL_COLUMN> void Recv_Columns(T_ARRAYS_HORIZONTAL_COLUMN& columns,const ARRAY<T_BOX_HORIZONTAL_INT>& regions,const int tag,
        const MPI::Status& probe_status) const;
//#####################################################################
};
}
#endif
#endif
