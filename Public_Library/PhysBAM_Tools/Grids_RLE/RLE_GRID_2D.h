//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_2D__
#define __RLE_GRID_2D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_RUN_2D.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
namespace PhysBAM{

template<class TV> class RANGE;
template<class T> class RLE_GRID_ITERATOR_CELL_2D;
template<class T> class RLE_GRID_ITERATOR_BLOCK_2D;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_2D_HELPER;
template<class T> class RLE_GRID_2D;
template<class T_GRID,int axis> class RLE_GRID_ITERATOR_FACE_HORIZONTAL;
template<class T_GRID> class RLE_GRID_ITERATOR_FACE_Y;
template<class TV> struct RLE_TAG;

template<class T>
class RLE_POLICY_2D
{
public:
    static const int dimension=2;
    static const int number_of_faces_per_cell=4;
    static const int number_of_neighbors_per_cell=4;
    static const int number_of_cells_per_block=4;
    static const int number_of_faces_per_block=12;

    typedef T SCALAR;
    typedef int INDEX;
    typedef VECTOR<T,dimension> VECTOR_T;
    typedef VECTOR<int,dimension> VECTOR_INT;
    typedef RLE_TAG<VECTOR_T> GRID_TAG;
    typedef RLE_GRID_ITERATOR_FACE_HORIZONTAL<RLE_GRID_2D<T>,1> FACE_X_ITERATOR;
    typedef RLE_GRID_ITERATOR_FACE_HORIZONTAL<RLE_GRID_2D<T>,3> FACE_Z_ITERATOR;
    typedef RLE_GRID_ITERATOR_FACE_Y<RLE_GRID_2D<T> > FACE_Y_ITERATOR;

    typedef RLE_RUN_2D RUN;
    typedef RLE_BLOCK_RUN_2D BLOCK_RUN;
    typedef RANGE<VECTOR<T,1> > BOX_HORIZONTAL;
    typedef RANGE<VECTOR<int,1> > BOX_HORIZONTAL_INT;
    typedef GRID<VECTOR<T,1> > HORIZONTAL_GRID;
    typedef VECTOR<T,1> VECTOR_HORIZONTAL;
    typedef VECTOR<int,1> VECTOR_HORIZONTAL_INT;
    typedef ARRAY<T,VECTOR<int,1> > ARRAYS_HORIZONTAL;
    typedef RLE_GRID_ITERATOR_BLOCK_2D<T> BLOCK;
    typedef RLE_GRID_ITERATOR_CELL_2D<T> CELL_ITERATOR;
    typedef RLE_GRID_ITERATOR_BLOCK_2D<T> BLOCK_ITERATOR;
    typedef LINEAR_INTERPOLATION_MAC_2D_HELPER<RLE_GRID_2D<T> > LINEAR_INTERPOLATION_MAC_HELPER;
    typedef HORIZONTAL_GRID PARALLEL_GRID;
};

template<class T>
class RLE_GRID_2D:public RLE_GRID<T,RLE_POLICY_2D<T> >
{
    typedef VECTOR<T,2> TV;
public:
    typedef RLE_GRID<T,RLE_POLICY_2D<T> > BASE;
    typedef typename BASE::CELL_ITERATOR CELL_ITERATOR;typedef typename BASE::FACE_X_ITERATOR FACE_X_ITERATOR;typedef typename BASE::FACE_Y_ITERATOR FACE_Y_ITERATOR;
    typedef typename BASE::FACE_Z_ITERATOR FACE_Z_ITERATOR;typedef typename BASE::BLOCK_ITERATOR BLOCK_ITERATOR;typedef typename BASE::RUN RUN;typedef typename BASE::BLOCK_RUN BLOCK_RUN;

    using BASE::uniform_grid;using BASE::number_of_ghost_cells;using BASE::minimum_vertical_space;using BASE::long_run_cells;using BASE::minimum_long_run_length; using BASE::negative_bandwidth;
    using BASE::positive_bandwidth;using BASE::domain;using BASE::columns;using BASE::jmin;using BASE::jmax;using BASE::horizontal_grid;using BASE::mpi_grid;using BASE::number_of_cells;
    using BASE::number_of_faces;using BASE::short_cell_neighbors;using BASE::short_face_neighbors;using BASE::Clean_Memory;using BASE::Clamped_Run_In_Column;using BASE::short_face_locations;
    using BASE::number_of_blocks;BASE::block_run_offsets;using BASE::block_runs;using BASE::last_face_numbers;using BASE::Block_Range_In_Column;

    RLE_GRID_2D();
    ~RLE_GRID_2D();

    const RUN* Clamped_Run_In_Column(const int i,const int j) const
    {return Clamped_Run_In_Column(columns(i),j);}

    void Clamped_Run(const VECTOR<T,2>& X,int& i,int& j,const RUN*& run) const
    {i=Clamped_i(X.x);j=Clamped_j(X.y);run=Clamped_Run_In_Column(i,j);}

    void Clamped_Block(const VECTOR<T,2>& X,int& j,const BLOCK_RUN*& run) const
    {int i=min(columns.domain.max_corner.x-2,columns.domain.min_corner.x+1+max(0,(int)((X.x-domain.min_corner.x)*uniform_grid.one_over_dX.x-(T).5-columns.domain.min_corner.x)));
    j=min(jmax-2,jmin+max(0,(int)((X.y-domain.min_corner.y)*uniform_grid.one_over_dX.y-(T).5)));
    run=Block_In_Column(j,block_run_offsets(i),block_run_offsets(i+1)-1);}    

    const BLOCK_RUN* Clamped_Block(const VECTOR<int,2>& I) const
    {return Block_In_Column(I.y,block_run_offsets(I.x),block_run_offsets(I.x+1)-1);}

    int Clamped_Block_Index(const VECTOR<T,2>& X,const int number_of_ghost_cells) const
    {int j;const RLE_BLOCK_RUN_2D *run;Clamped_Block(X,j,run);
    if(run && run->i>=1-number_of_ghost_cells && run->i<=uniform_grid.counts.x-2+number_of_ghost_cells) return run->block+j-run->jmin;
    return 0;}

    RANGE<VECTOR<int,1> > Block_Range(const VECTOR<int,1>& I,const int jmin,const int jmax) const
    {return Block_Range_In_Column(jmin,jmax,block_run_offsets(I.x),block_run_offsets(I.x+1)-1);}

    template<class T2> T2 Cell_Value(const ARRAY<T2>& u,const VECTOR<int,2>& I) const
    {return Cell_Value(u,I.x,I.y);}

    template<class T2> T2 Cell_Value(const ARRAY<T2>& u,const int i,const int j) const
    {assert(long_run_cells==2);const RLE_RUN_2D* run=Clamped_Run_In_Column(i,j);if(!run->is_long) return u(run->cell+j-run->jmin);
    return LINEAR_INTERPOLATION<T,T2>::Linear(u(run->cell),u(run->cell+1),(T)(j-run->jmin)/((run+1)->jmin-1-run->jmin));}

    int Face_Axis(const int face) const
    {return face<=last_face_numbers.y?2:1;}

    RANGE<VECTOR<int,1> > Cells_In_Slice(const int i) const
    {return RANGE<VECTOR<int,1> >(columns(i)(1).cell,columns(i+1)(1).cell-1);}

    RANGE<VECTOR<int,1> > Cells_In_Slices(const int i_start,const int i_end) const
    {return RANGE<VECTOR<int,1> >(columns(i_start)(1).cell,columns(i_end+1)(1).cell-1);}

    template<class F,class T1> static void Horizontal_Face_Loop(T1& a1)
    {F::template Apply<FACE_X_ITERATOR>(a1);}

    template<class F,class T1,class T2> static void Horizontal_Face_Loop(T1& a1,T2& a2)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2);}

    template<class F,class T1,class T2,class T3> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3);}

    template<class F,class T1,class T2,class T3,class T4> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4);}

    template<class F,class T1,class T2,class T3,class T4,class T5> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4,a5);}

    template<class F,class T1,class T2,class T3,class T4,class T5,class T6> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5,T6& a6)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4,a5,a6);}

    template<class F,class T1,class T2,class T3,class T4,class T5,class T6,class T7> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5,T6& a6,T7& a7)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4,a5,a6,a7);}

//#####################################################################
    void Initialize(const RLE_INITIAL_PHI_HELPER<T,2>& initial_phi,const ARRAY<int,VECTOR<int,1> >* ground_j);
    void Adjust_Vertical_Space();
    static void Rebuild(const RLE_GRID_2D<T>& old_grid,RLE_GRID_2D<T>& new_grid,const ARRAY<bool>& cell_should_be_long,const int extra_cells,const ARRAY<int,VECTOR<int,1> >* ground_j);
    void Compute_Auxiliary_Information();
    void Find_Ghost_Regions(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& sentinels,const int number_of_ghost_cells,const bool include_corners=true) const;
private:
    void Compute_Block_Columns();
    void Calculate_Short_Cell_Neighbors() const;
    void Calculate_Short_Face_Neighbors() const;
    void Calculate_Short_Face_Locations() const;
    friend class RLE_GRID<T,RLE_POLICY_2D<T> >;
//#####################################################################
};
}
#endif
#endif
