//#####################################################################
// Copyright 2005, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_CELL_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATOR_CELL_2D__
#define __RLE_GRID_ITERATOR_CELL_2D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_IN_COLUMN.h>
namespace PhysBAM{

template<class T_input>
class RLE_GRID_ITERATOR_CELL_2D:public RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T_input> >
{
    typedef T_input T;
public:
    typedef RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> > BASE;
    using BASE::grid;using BASE::column;using BASE::j;using BASE::Short;using BASE::jmax;using BASE::Initialize_Column;

    int i_end,i;

    RLE_GRID_ITERATOR_CELL_2D(const RLE_GRID_2D<T>& grid_input,const int number_of_ghost_cells,const RANGE<VECTOR<int,1> >& sentinels=(RANGE<VECTOR<int,1> >(0,0)),const bool top_sentinels=false)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> >(grid_input,top_sentinels)
    {
        Initialize(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices().Thickened(number_of_ghost_cells)+sentinels);
    }

    RLE_GRID_ITERATOR_CELL_2D(const RLE_GRID_2D<T>& grid_input,const int number_of_ghost_cells,const RANGE<VECTOR<int,1> >& sentinels,const int bottom_sentinels,const int top_sentinels)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> >(grid_input,bottom_sentinels,top_sentinels)
    {
        Initialize(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices().Thickened(number_of_ghost_cells)+sentinels);
    }

    RLE_GRID_ITERATOR_CELL_2D(const RLE_GRID_2D<T>& grid_input,const RANGE<VECTOR<int,1> >& region,const bool top_sentinels=false)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> >(grid_input,top_sentinels)
    {
        Initialize(region);
    }

    RLE_GRID_ITERATOR_CELL_2D(const RLE_GRID_2D<T>& grid_input,const RANGE<VECTOR<int,1> >& region,const int bottom_sentinels,const int top_sentinels)
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> >(grid_input,bottom_sentinels,top_sentinels)
    {
        Initialize(region);
    }

    RLE_GRID_ITERATOR_CELL_2D(const RLE_GRID_2D<T>& grid_input,const VECTOR<T,2>& X) // creates an iterator that can't be incremented
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> >(grid_input,false)
    {
        int i_input,j_input;const RLE_RUN_2D* run_input;grid.Clamped_Run(X,i_input,j_input,run_input);
        Initialize(i_input,j_input,run_input);
    }

    RLE_GRID_ITERATOR_CELL_2D(const RLE_GRID_2D<T>& grid_input,const VECTOR<int,2>& I) // creates an iterator that can't be incremented
        :RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> >(grid_input,false)
    {
        Initialize(I.x,I.y,grid.Clamped_Run_In_Column(I));
    }

private:
    void Initialize(const RANGE<VECTOR<int,1> >& region)
    {i=region.min_corner.x;i_end=region.max_corner.x;column=&grid.columns(i);Initialize_Column();}

    void Initialize(const int i_input,const int j_input,const RLE_RUN_2D* run_input)
    {i_end=i=i_input;Initialize_Column(j_input,run_input);}
public:

    operator bool() const
    {return i<=i_end;}

    void operator++(int)
    {assert(*this);
    if(!RLE_GRID_ITERATOR_CELL_IN_COLUMN<RLE_GRID_2D<T> >::operator++()){
        if(++i>i_end) return; // finished
        column=&grid.columns(i);Initialize_Column();}} // start a new column

    T x() const
    {return grid.uniform_grid.Axis_X_plus_half(i,1);}

    VECTOR<int,2> I() const
    {return VECTOR<int,2>(i,j);}

    VECTOR<int,1> Horizontal_Index() const
    {return VECTOR<int,1>(i);}

    VECTOR<T,2> X() const
    {assert(Short());return grid.uniform_grid.Center(i,j);}

    VECTOR<T,2> Face_Location(const int face) const
    {assert(Short() && 1<=face&&face<=4);
    VECTOR<T,2> face_X=X();int axis=(face-1)/2+1;T dx_over_two=(T).5*grid.uniform_grid.dX[axis];
    face_X[axis]+=face&1?-dx_over_two:dx_over_two;return face_X;}

    VECTOR<T,2> Center() const
    {return Short()?X():VECTOR<T,2>(grid.uniform_grid.Axis_X_plus_half(i,1),(T).5*(grid.uniform_grid.Axis_X(j,2)+grid.uniform_grid.Axis_X(jmax(),2)));}

    static RANGE<VECTOR<int,1> > Indices_In_Slices(const RLE_GRID_2D<T>& grid,const int i_start,const int i_end)
    {return grid.Cells_In_Slices(i_start,i_end);}

//#####################################################################
};
}
#endif
#endif
