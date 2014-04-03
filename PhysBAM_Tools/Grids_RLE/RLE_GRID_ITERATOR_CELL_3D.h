//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_CELL_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATORS_3D__
#define __RLE_GRID_ITERATORS_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_IN_COLUMN.h>
namespace PhysBAM{

template<class T_input>
class RLE_GRID_ITERATOR_CELL_3D:public RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T_input> >
{
    typedef T_input T;
public:
    typedef RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> > BASE;
    using BASE::grid;using BASE::column;using BASE::j;using BASE::Short;using BASE::jmax;using BASE::Initialize_Column;

    int i_end,i,ij_start,ij_end,ij;

    RLE_GRID_ITERATOR_CELL_3D(const RLE_GRID_3D<T>& grid_input,const int number_of_ghost_cells,const RANGE<VECTOR<int,2> >& sentinels=(RANGE<VECTOR<int,2> >(0,0,0,0)),const bool top_sentinels=false)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> >(grid_input,top_sentinels)
    {
        Initialize(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices().Thickened(number_of_ghost_cells)+sentinels);
    }

    RLE_GRID_ITERATOR_CELL_3D(const RLE_GRID_3D<T>& grid_input,const int number_of_ghost_cells,const RANGE<VECTOR<int,2> >& sentinels,const int bottom_sentinels,const int top_sentinels)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> >(grid_input,bottom_sentinels,top_sentinels)
    {
        Initialize(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices().Thickened(number_of_ghost_cells)+sentinels);
    }

    RLE_GRID_ITERATOR_CELL_3D(const RLE_GRID_3D<T>& grid_input,const RANGE<VECTOR<int,2> >& region,const bool top_sentinels=false)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> >(grid_input,top_sentinels)
    {
        Initialize(region);
    }

    RLE_GRID_ITERATOR_CELL_3D(const RLE_GRID_3D<T>& grid_input,const RANGE<VECTOR<int,2> >& region,const int bottom_sentinels,const int top_sentinels)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> >(grid_input,bottom_sentinels,top_sentinels)
    {
        Initialize(region);
    }

    RLE_GRID_ITERATOR_CELL_3D(const RLE_GRID_3D<T>& grid_input,const VECTOR<T,3>& X) // creates an iterator that can't be incremented
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> >(grid_input,false)
    {
        int i_input,j_input,ij_input;const RLE_RUN_3D* run_input;grid.Clamped_Run(X,i_input,j_input,ij_input,run_input);
        Initialize(i_input,j_input,ij_input,run_input);
    }

    RLE_GRID_ITERATOR_CELL_3D(const RLE_GRID_3D<T>& grid_input,const VECTOR<int,3>& I) // creates an iterator that can't be incremented
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> >(grid_input,false)
    {
        Initialize(I.x,I.y,I.z,grid.Clamped_Run_In_Column(I));
    }

private:
    void Initialize(const RANGE<VECTOR<int,2> >& region)
    {i=region.min_corner.x;i_end=region.max_corner.x;ij_start=ij=region.min_corner.y;ij_end=region.max_corner.y;column=&grid.columns(i,ij);Initialize_Column();}

    void Initialize(const int i_input,const int j_input,const int ij_input,const RLE_RUN_3D* run_input)
    {i_end=i=i_input;ij_start=ij_end=ij=ij_input;Initialize_Column(j_input,run_input);}
public:

    operator bool() const
    {return i<=i_end;}

    void operator++(int)
    {assert(*this);
    if(!RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_3D<T> >::operator++()){
        if(++ij>ij_end){
            ij=ij_start;
            if(++i>i_end) return;} // finished
        column=&grid.columns(i,ij);Initialize_Column();}} // start a new column

    T x() const
    {return grid.uniform_grid.Axis_X_plus_half(i,1);}

    T z() const
    {return grid.uniform_grid.Axis_X_plus_half(ij,3);}

    VECTOR<int,3> I() const
    {return VECTOR<int,3>(i,j,ij);}

    VECTOR<int,2> Horizontal_Index() const
    {return VECTOR<int,2>(i,ij);}

    VECTOR<T,3> X() const
    {assert(Short());return grid.uniform_grid.Center(i,j,ij);}

    VECTOR<T,3> Face_Location(const int face) const
    {assert(Short() && 1<=face&&face<=6);
    VECTOR<T,3> face_X=X();int axis=(face-1)/2+1;T dx_over_two=(T).5*grid.uniform_grid.dX[axis];
    face_X[axis]+=face&1?-dx_over_two:dx_over_two;return face_X;}

    VECTOR<T,3> Center() const
    {return Short()?X():VECTOR<T,3>(grid.uniform_grid.Axis_X_plus_half(i,1),(T).5*(grid.uniform_grid.Axis_X(j,2)+grid.uniform_grid.Axis_X(jmax(),2)),grid.uniform_grid.Axis_X_plus_half(ij,3));}

    static RANGE<VECTOR<int,1> > Indices_In_Slice(const RLE_GRID_3D<T>& grid,const int i,const int ij_start,const int ij_end)
    {return grid.Cells_In_Slice(i,ij_start,ij_end);}

    static RANGE<VECTOR<int,1> > Indices_In_Slice(const RLE_GRID_3D<T>& grid,const int i)
    {return grid.Cells_In_Slice(i);}

//#####################################################################
};
}
#endif
#endif
