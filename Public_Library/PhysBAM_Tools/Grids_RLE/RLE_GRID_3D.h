//#####################################################################
// Copyright 2005, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_3D__
#define __RLE_GRID_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_RUN_3D.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
namespace PhysBAM{

template<class T> class RLE_GRID_ITERATOR_CELL_3D;
template<class T> class RLE_GRID_ITERATOR_BLOCK_3D;
template<class T> class RLE_GRID_TRANSFER_ITERATOR_BLOCK_3D;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_3D_HELPER;
template<class T> class RLE_GRID_3D;
template<class T_GRID,int axis> class RLE_GRID_ITERATOR_FACE_HORIZONTAL;
template<class T_GRID> class RLE_GRID_ITERATOR_FACE_Y;
template<class TV> struct RLE_TAG;

template<class T>
class RLE_POLICY_3D
{
public:
    static const int dimension=3;
    static const int number_of_faces_per_cell=6;
    static const int number_of_neighbors_per_cell=6;
    static const int number_of_cells_per_block=8;
    static const int number_of_faces_per_block=36;

    typedef T SCALAR;
    typedef int INDEX;
    typedef VECTOR<T,dimension> VECTOR_T;
    typedef VECTOR<int,dimension> VECTOR_INT;
    typedef RLE_TAG<VECTOR_T> GRID_TAG;
    typedef RLE_GRID_ITERATOR_FACE_HORIZONTAL<RLE_GRID_3D<T>,1> FACE_X_ITERATOR;
    typedef RLE_GRID_ITERATOR_FACE_HORIZONTAL<RLE_GRID_3D<T>,3> FACE_Z_ITERATOR;
    typedef RLE_GRID_ITERATOR_FACE_Y<RLE_GRID_3D<T> > FACE_Y_ITERATOR;

    typedef RLE_RUN_3D RUN;
    typedef RLE_BLOCK_RUN_3D BLOCK_RUN;
    typedef RANGE<VECTOR<T,2> > BOX_HORIZONTAL;
    typedef RANGE<VECTOR<int,2> > BOX_HORIZONTAL_INT;
    typedef GRID<VECTOR<T,2> > HORIZONTAL_GRID;
    typedef VECTOR<T,2> VECTOR_HORIZONTAL;
    typedef VECTOR<int,2> VECTOR_HORIZONTAL_INT;
    typedef ARRAY<T,VECTOR<int,2> > ARRAYS_HORIZONTAL;
    typedef RLE_GRID_ITERATOR_BLOCK_3D<T> BLOCK;
    typedef RLE_GRID_ITERATOR_CELL_3D<T> CELL_ITERATOR;
    typedef RLE_GRID_ITERATOR_BLOCK_3D<T> BLOCK_ITERATOR;
    typedef RLE_GRID_TRANSFER_ITERATOR_BLOCK_3D<T> BLOCK_TRANSFER_ITERATOR;
    typedef LINEAR_INTERPOLATION_MAC_3D_HELPER<RLE_GRID_3D<T> > LINEAR_INTERPOLATION_MAC_HELPER;
    typedef HORIZONTAL_GRID PARALLEL_GRID;
};

template<class T>
class RLE_GRID_3D:public RLE_GRID<T,RLE_POLICY_3D<T> >
{
    typedef VECTOR<T,3> TV;
public:
    typedef RLE_GRID<T,RLE_POLICY_3D<T> > BASE;
    typedef typename BASE::BLOCK_RUN BLOCK_RUN;typedef typename BASE::CELL_ITERATOR CELL_ITERATOR;typedef typename BASE::FACE_X_ITERATOR FACE_X_ITERATOR;
    typedef typename BASE::FACE_Y_ITERATOR FACE_Y_ITERATOR;typedef typename BASE::FACE_Z_ITERATOR FACE_Z_ITERATOR;typedef typename BASE::BLOCK_ITERATOR BLOCK_ITERATOR;

    using BASE::uniform_grid;using BASE::number_of_ghost_cells;using BASE::minimum_vertical_space;using BASE::minimum_long_run_length;using BASE::long_run_cells;using BASE::negative_bandwidth;
    using BASE::positive_bandwidth;using BASE::number_of_cells;using BASE::number_of_faces;using BASE::horizontal_grid;using BASE::domain;using BASE::columns;using BASE::mpi_grid;
    using BASE::jmin;using BASE::jmax;using BASE::short_cell_neighbors;using BASE::short_face_neighbors;using BASE::short_face_locations;using BASE::Face_Loop;using BASE::last_face_numbers;
    using BASE::Clean_Memory;using BASE::Clamped_Run_In_Column;using BASE::number_of_blocks;using BASE::block_run_offsets;using BASE::block_runs;typedef typename BASE::RUN RUN;
    using BASE::Block_Range_In_Column;

    RLE_GRID_3D();
    ~RLE_GRID_3D();

    int Clamped_ij(const T z) const
    {return min(columns.domain.max_corner.y-1,columns.domain.min_corner.y+1+max(0,(int)((z-domain.min_corner.z)*uniform_grid.one_over_dX.z-columns.domain.min_corner.y)));}

    const RUN* Clamped_Run_In_Column(const int i,const int j,const int ij) const
    {return Clamped_Run_In_Column(columns(i,ij),j);}

    void Clamped_Run(const VECTOR<T,3>& X,int& i,int& j,int& ij,const RUN*& run) const
    {i=Clamped_i(X.x);j=Clamped_j(X.y);ij=Clamped_ij(X.z);run=Clamped_Run_In_Column(i,j,ij);}

    int Block_i(const T x) const
    {return min(columns.domain.max_corner.x-2,columns.domain.min_corner.x+1+max(0,(int)((x-domain.min_corner.x)*uniform_grid.one_over_dX.x-(T).5-columns.domain.min_corner.x)));}

    int Block_ij(const T z) const
    {return min(columns.domain.max_corner.y-2,columns.domain.min_corner.y+1+max(0,(int)((z-domain.min_corner.z)*uniform_grid.one_over_dX.z-(T).5-columns.domain.min_corner.y)));}

    void Clamped_Block(const VECTOR<T,3>& X,int& i,int& j,const BLOCK_RUN*& run) const
    {i=Block_i(X.x);int ij=Block_ij(X.z);
    j=min(jmax-2,jmin+max(0,(int)((X.y-domain.min_corner.y)*uniform_grid.one_over_dX.y-(T).5)));
    run=Block_In_Column(j,block_run_offsets(i,ij),block_run_offsets(i,ij+1)-1);}

    const BLOCK_RUN* Clamped_Block(const VECTOR<int,3>& I) const
    {return Block_In_Column(I.y,block_run_offsets(I.x,I.z),block_run_offsets(I.x,I.z+1)-1);}

    int Clamped_Block_Index(const VECTOR<T,3>& X,const int number_of_ghost_cells) const
    {int i,j;const BLOCK_RUN *run;Clamped_Block(X,i,j,run);
    if(run && i>=1-number_of_ghost_cells && i<=uniform_grid.counts.x-2+number_of_ghost_cells && run->ij>=1-number_of_ghost_cells && run->ij<=uniform_grid.counts.z-2+number_of_ghost_cells)
        return run->block+j-run->jmin;
    return 0;}

    RANGE<VECTOR<int,1> > Block_Range(const VECTOR<int,2>& I,const int jmin,const int jmax) const
    {return Block_Range_In_Column(jmin,jmax,block_run_offsets(I.x,I.y),block_run_offsets(I.x,I.y+1)-1);}

    template<class T2> T2 Cell_Value(const ARRAY<T2>& u,const VECTOR<int,3>& I) const
    {return Cell_Value(u,I.x,I.y,I.z);}

    template<class T2> T2 Cell_Value(const ARRAY<T2>& u,const int i,const int j,const int ij) const
    {assert(long_run_cells==2);const RUN* run=Clamped_Run_In_Column(i,j,ij);if(!run->is_long) return u(run->cell+j-run->jmin);
    return LINEAR_INTERPOLATION<T,T2>::Linear(u(run->cell),u(run->cell+1),(T)(j-run->jmin)/((run+1)->jmin-1-run->jmin));}

    int Face_Axis(const int face) const
    {return face<=last_face_numbers.y?2:face<=last_face_numbers.x?1:3;}

    RANGE<VECTOR<int,1> > Cells_In_Slice(const int i,const int ij_start,const int ij_end) const
    {return RANGE<VECTOR<int,1> >(columns(i,ij_start)(1).cell,columns(i,ij_end+1)(1).cell-1);}

    RANGE<VECTOR<int,1> > Cells_In_Slice(const int i) const
    {return Cells_In_Slice(i,columns.domain.min_corner.y+1,columns.domain.max_corner.y-1);}

    template<class F,class T1> static void Horizontal_Face_Loop(T1& a1)
    {F::template Apply<FACE_X_ITERATOR>(a1);F::template Apply<FACE_Z_ITERATOR>(a1);}

    template<class F,class T1,class T2> static void Horizontal_Face_Loop(T1& a1,T2& a2)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2);F::template Apply<FACE_Z_ITERATOR>(a1,a2);}

    template<class F,class T1,class T2,class T3> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3);F::template Apply<FACE_Z_ITERATOR>(a1,a2,a3);}

    template<class F,class T1,class T2,class T3,class T4> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4);F::template Apply<FACE_Z_ITERATOR>(a1,a2,a3,a4);}

    template<class F,class T1,class T2,class T3,class T4,class T5> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4,a5);F::template Apply<FACE_Z_ITERATOR>(a1,a2,a3,a4,a5);}

    template<class F,class T1,class T2,class T3,class T4,class T5,class T6> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5,T6& a6)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4,a5,a6);F::template Apply<FACE_Z_ITERATOR>(a1,a2,a3,a4,a5,a6);}

    template<class F,class T1,class T2,class T3,class T4,class T5,class T6,class T7> static void Horizontal_Face_Loop(T1& a1,T2& a2,T3& a3,T4& a4,T5& a5,T6& a6,T7& a7)
    {F::template Apply<FACE_X_ITERATOR>(a1,a2,a3,a4,a5,a6,a7);F::template Apply<FACE_Z_ITERATOR>(a1,a2,a3,a4,a5,a6,a7);}

//#####################################################################
    void Initialize(const RLE_INITIAL_PHI_HELPER<T,3>& initial_phi,const ARRAY<int,VECTOR<int,2> >* ground_j,const bool verbose=false);
    void Adjust_Vertical_Space();
    void Adjust_Vertical_Space(const int max_jmin,const int min_jmax);
    static void Rebuild(const RLE_GRID_3D<T>& old_grid,RLE_GRID_3D<T>& new_grid,const ARRAY<bool>& cell_should_be_long,const int extra_cells,const ARRAY<int,VECTOR<int,2> >* ground_j);
    void Compute_Auxiliary_Information();
    void Find_Ghost_Regions(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& sentinels,const int number_of_ghost_cells,const bool include_corners=true,
        const bool use_separate_regions_for_corners=false) const;
private:
    void Compute_Block_Columns();
    void Calculate_Short_Cell_Neighbors() const;
    void Calculate_Short_Face_Neighbors() const;
    void Calculate_Short_Face_Locations() const;
    friend class RLE_GRID<T,RLE_POLICY_3D<T> >;
//#####################################################################
};
}
#endif
#endif
