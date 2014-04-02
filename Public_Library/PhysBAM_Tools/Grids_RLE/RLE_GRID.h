//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID__
#define __RLE_GRID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_BLOCK_RUNS.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T_GRID,class T_ITERATOR> class RLE_GRID_SIMPLE_ITERATOR;
template<class T_GRID> class MPI_RLE_GRID;
    
template<class T, int d>
class RLE_INITIAL_PHI_HELPER
{
public:
    typedef VECTOR<T,d> TV;
    virtual ~RLE_INITIAL_PHI_HELPER() {};
    virtual T operator()(const TV& v) const=0;
};

template<class T,class POLICY>
class RLE_GRID:public POLICY,public NONCOPYABLE
{
    typedef typename POLICY::VECTOR_T TV;typedef typename POLICY::VECTOR_INT TV_INT;typedef typename RLE_GRID_POLICY<TV>::RLE_GRID T_GRID;
    typedef typename POLICY::RUN T_RUN;typedef typename POLICY::BLOCK_RUN T_BLOCK_RUN;
public:
    typedef typename POLICY::HORIZONTAL_GRID T_HORIZONTAL_GRID;typedef typename POLICY::CELL_ITERATOR CELL_ITERATOR;typedef typename POLICY::FACE_Y_ITERATOR FACE_Y_ITERATOR;
    typedef typename POLICY::BLOCK_ITERATOR BLOCK_ITERATOR;typedef typename POLICY::ARRAYS_HORIZONTAL T_ARRAYS_HORIZONTAL;typedef MPI_RLE_GRID<T_GRID> T_MPI_GRID;
    typedef typename REBIND<T_ARRAYS_HORIZONTAL,ARRAY<T_RUN> >::TYPE T_ARRAYS_HORIZONTAL_ARRAY_RUN;typedef typename POLICY::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename REBIND<T_ARRAYS_HORIZONTAL,int>::TYPE T_ARRAYS_HORIZONTAL_INT;typedef typename T_HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;
    typedef typename POLICY::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;

    GRID<TV> uniform_grid;
    int number_of_ghost_cells;
    int minimum_vertical_space,minimum_long_run_length;
    int long_run_cells,long_run_faces_horizontal;
    T negative_bandwidth,positive_bandwidth;
    int number_of_cells,number_of_faces,number_of_blocks;
    T_HORIZONTAL_GRID horizontal_grid;
    RANGE<TV> domain;
    T_ARRAYS_HORIZONTAL_ARRAY_RUN columns; // last run in each column is a long sentinel, as are first and last columns
    T_ARRAYS_HORIZONTAL_INT block_run_offsets;
    ARRAY<T_BLOCK_RUN> block_runs;
    T_MPI_GRID* mpi_grid;
    int jmin,jmax;
protected:
    TV_INT last_face_numbers;
    mutable ARRAY<VECTOR<int,POLICY::number_of_neighbors_per_cell> > short_cell_neighbors;
    mutable ARRAY<VECTOR<int,POLICY::number_of_neighbors_per_cell> > short_face_neighbors;
    mutable ARRAY<TV> short_face_locations;

    RLE_GRID();
    ~RLE_GRID();
public:

    void Set_Uniform_Grid(const GRID<TV>& uniform_grid_input,const int number_of_ghost_cells_input=3,const int minimum_vertical_space_input=10)
    {if(columns.counts.x) PHYSBAM_FATAL_ERROR();
    if(uniform_grid_input.Is_MAC_Grid() || number_of_ghost_cells_input<0 || minimum_vertical_space_input<0) PHYSBAM_FATAL_ERROR();
    if(!uniform_grid_input.Is_Isotropic()) PHYSBAM_FATAL_ERROR();
    uniform_grid=uniform_grid_input;number_of_ghost_cells=number_of_ghost_cells_input;minimum_vertical_space=minimum_vertical_space_input;
    assert(Minimum_Horizontal_Cells(uniform_grid)>number_of_ghost_cells);
    horizontal_grid=uniform_grid.Get_Horizontal_Grid();}

    static int Minimum_Horizontal_Cells(const GRID<TV>& uniform_grid)
    {return uniform_grid.numbers_of_cells.Min();}

    void Set_Minimum_Long_Run_Length(const int minimum_long_run_length_input=2)
    {minimum_long_run_length=minimum_long_run_length_input;assert(1<=minimum_long_run_length && 2*minimum_long_run_length<minimum_vertical_space);}

    void Set_Long_Run_Sample_Counts(const int long_run_cells_input=2,const int long_run_faces_horizontal_input=1)
    {long_run_cells=long_run_cells_input;long_run_faces_horizontal=long_run_faces_horizontal_input;
    assert(1<=long_run_cells && long_run_cells<=2 && 1<=long_run_faces_horizontal && long_run_faces_horizontal<=2);}

    void Set_Linear_Pressure_And_Constant_Horizontal_Velocity_Profile()
    {Set_Minimum_Long_Run_Length(2);Set_Long_Run_Sample_Counts(2,1);}

    void Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile()
    {Set_Minimum_Long_Run_Length(2);Set_Long_Run_Sample_Counts(2,2);}

    void Set_Positive_Bandwidth_In_Cells(const T positive_bandwidth_in_cells=(T)3)
    {T max_length=uniform_grid.dX.Max();if(!max_length || positive_bandwidth_in_cells<(T)2) PHYSBAM_FATAL_ERROR();
    positive_bandwidth=positive_bandwidth_in_cells*max_length;}

    void Set_Negative_Bandwidth_In_Cells(const T negative_bandwidth_in_cells=(T)5)
    {T max_length=uniform_grid.dX.Max();if(!max_length || negative_bandwidth_in_cells<(T)2) PHYSBAM_FATAL_ERROR();
    negative_bandwidth=negative_bandwidth_in_cells*max_length;}

    void Set_Negative_Bandwidth_From_Optical_Depth(const T optical_depth)
    {T max_length=uniform_grid.dX.Max();if(!max_length || optical_depth<(T)3*max_length) PHYSBAM_FATAL_ERROR();
    negative_bandwidth=optical_depth;}

    T Minimum_Edge_Length() const
    {return uniform_grid.Minimum_Edge_Length();}

    const RANGE<TV>& Domain() const
    {return domain;}

    TV Clamp(const TV& X) const
    {return domain.Clamp(X);}

    bool Outside(const TV& X) const
    {return domain.Lazy_Outside(X);}

    int Clamped_i(const T x) const
    {return min(columns.domain.max_corner.x-1,columns.domain.min_corner.x+1+max(0,(int)((x-domain.min_corner.x)*uniform_grid.one_over_dX.x-columns.domain.min_corner.x)));}

    int Clamped_j(const T y) const
    {return min(jmax-1,jmin+max(0,(int)((y-domain.min_corner.y)*uniform_grid.one_over_dX.y)));}

    const T_RUN* Clamped_Run_In_Column(const TV_INT& I) const
    {return Clamped_Run_In_Column(columns(I.Horizontal_Vector()),I.y);}

protected:
    const T_RUN* Clamped_Run_In_Column(const ARRAY<T_RUN>& column,const int j) const
    {int lo=1,hi=column.m-1;
    while(lo<hi){int mid=(lo+hi+1)>>1;
        if(j<column(mid).jmin) hi=mid-1;else lo=mid;}
    return &column(lo);}

    int Block_Run_Index_In_Column(const int j,const int block_lo,const int block_hi) const
    {int lo=block_lo,hi=block_hi;
    while(lo<hi){int mid=(lo+hi+1)>>1;
        if(j<block_runs(mid).jmin) hi=mid-1;else lo=mid;}
    return lo;}

    const T_BLOCK_RUN* Block_In_Column(const int j,const int block_lo,const int block_hi) const
    {if(block_hi<block_lo) return 0;
    int index=Block_Run_Index_In_Column(j,block_lo,block_hi);
    if(block_runs(index).jmin<=j && j<=block_runs(index).jmax) return &block_runs(index);else return 0;}

    RANGE<VECTOR<int,1> > Block_Range_In_Column(const int jmin,const int jmax,const int block_lo,const int block_hi) const
    {if(block_hi<block_lo) return RANGE<VECTOR<int,1> >(1,0);
    int hi_run=Block_Run_Index_In_Column(jmax,block_lo,block_hi);
    if(jmax<block_runs(hi_run).jmin) return RANGE<VECTOR<int,1> >(1,0);
    int lo_run=hi_run;while(block_lo<lo_run && jmin<=block_runs(lo_run-1).jmax) lo_run--;
    if(jmin>block_runs(lo_run).jmax) return RANGE<VECTOR<int,1> >(1,0);
    int first_block=block_runs(lo_run).block+max(0,jmin-block_runs(lo_run).jmin);
    int last_block=block_runs(hi_run).block+min((int)block_runs(hi_run).jmax,jmax)-block_runs(hi_run).jmin;
    return RANGE<VECTOR<int,1> >(first_block,last_block);}
public:

    const ARRAY<VECTOR<int,POLICY::number_of_neighbors_per_cell> >& Short_Cell_Neighbors() const
    {if(!short_cell_neighbors.m) static_cast<const T_GRID&>(*this).Calculate_Short_Cell_Neighbors();return short_cell_neighbors;}

    const ARRAY<VECTOR<int,POLICY::number_of_neighbors_per_cell> >& Short_Face_Neighbors() const
    {if(!short_face_neighbors.m) static_cast<const T_GRID&>(*this).Calculate_Short_Face_Neighbors();return short_face_neighbors;}

    const ARRAY<TV>& Short_Face_Locations() const
    {if(!short_face_locations.m) static_cast<const T_GRID&>(*this).Calculate_Short_Face_Locations();return short_face_locations;}

    template<class F,class T1> static void Face_Loop(T1& a1)
    {T_GRID::template Horizontal_Face_Loop<F>(a1);F::template Apply<FACE_Y_ITERATOR>(a1);}

    template<class F,class T1,class T2> static void Face_Loop(T1& a1,T2& a2)
    {T_GRID::template Horizontal_Face_Loop<F>(a1,a2);F::template Apply<FACE_Y_ITERATOR>(a1,a2);}

    template<class F,class T1,class T2,class T3> static void Face_Loop(T1& a1,T2& a2,T3& a3)
    {T_GRID::template Horizontal_Face_Loop<F>(a1,a2,a3);F::template Apply<FACE_Y_ITERATOR>(a1,a2,a3);}

    template<class F,class T1,class T2,class T3,class T4> static void Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4)
    {T_GRID::template Horizontal_Face_Loop<F>(a1,a2,a3,a4);F::template Apply<FACE_Y_ITERATOR>(a1,a2,a3,a4);}

    template<class F,class T1,class T2,class T3,class T4,class T5> static void Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5)
    {T_GRID::template Horizontal_Face_Loop<F>(a1,a2,a3,a4,a5);F::template Apply<FACE_Y_ITERATOR>(a1,a2,a3,a4,a5);}

    template<class F,class T1,class T2,class T3,class T4,class T5,class T6> static void Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5,T6& a6)
    {T_GRID::template Horizontal_Face_Loop<F>(a1,a2,a3,a4,a5,a6);F::template Apply<FACE_Y_ITERATOR>(a1,a2,a3,a4,a5,a6);}

    template<class F,class T1,class T2,class T3,class T4,class T5,class T6,class T7> static void Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5,T6& a6,T7& a7)
    {T_GRID::template Horizontal_Face_Loop<F>(a1,a2,a3,a4,a5,a6,a7);F::template Apply<FACE_Y_ITERATOR>(a1,a2,a3,a4,a5,a6,a7);}

    int Cell_Indices(const int ghost_cells=0) const
    {return number_of_cells;}

    int Block_Indices(const int ghost_cells=0) const
    {return number_of_blocks;}

    int Face_Indices(const int ghost_cells=0) const
    {return number_of_faces;}

    int Domain_Indices(const int ghost_cells=0) const
    {return number_of_cells;}

//#####################################################################
    void Clean_Memory();
    static void Transfer(RLE_GRID<T,POLICY>& old_grid,RLE_GRID<T,POLICY>& new_grid); // destroys old_grid
    void Topology_Changed();
    static void Transfer_Noncolumn_Data(const RLE_GRID<T,POLICY>& old_grid,RLE_GRID<T,POLICY>& new_grid);
    static void Union(const T_GRID& grid1,const T_GRID& grid2,T_GRID& union_grid,const T_ARRAYS_HORIZONTAL_INT& ground_j);
    template<class T2> void Put_Ghost_Cells(const T2& constant,ARRAY<T2>& array) const;
    template<class T2> void Put_Ghost_Faces(const T2& constant,ARRAY<T2>& array) const;
//#####################################################################
};
}
#endif
#endif
