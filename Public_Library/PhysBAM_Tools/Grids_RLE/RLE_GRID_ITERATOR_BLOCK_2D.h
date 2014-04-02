//#####################################################################
// Copyright 2005, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_BLOCK_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATOR_BLOCK_2D__
#define __RLE_GRID_ITERATOR_BLOCK_2D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK.h>
namespace PhysBAM{

template<class T_input>
class RLE_GRID_ITERATOR_BLOCK_2D:public RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T_input> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;
public:
    typedef RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T> > BASE;
    using BASE::grid;using BASE::run_end;using BASE::run;using BASE::j_end;using BASE::j;using BASE::dj;using BASE::Face_X;using BASE::Face_Y;

    RLE_GRID_ITERATOR_BLOCK_2D(const RLE_GRID_2D<T>& grid_input,const int number_of_ghost_cells)
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T> >(grid_input)
    {
        RANGE<VECTOR<int,1> > region=grid.horizontal_grid.Get_MAC_Grid().Domain_Indices()+Sentinels();
        Initialize(region.Thickened(number_of_ghost_cells-1));
    }

    RLE_GRID_ITERATOR_BLOCK_2D(const RLE_GRID_2D<T>& grid_input,const RANGE<VECTOR<int,1> >& region)
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T> >(grid_input)
    {
        Initialize(region);
    }

    RLE_GRID_ITERATOR_BLOCK_2D(const RLE_GRID_2D<T>& grid_input,const TV& X) // creates an iterator that cannot be incremented
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T> >(grid_input)
    {
        Initialize(X);
    }

    RLE_GRID_ITERATOR_BLOCK_2D(const RLE_GRID_2D<T>& grid_input,const VECTOR<int,2>& I) // creates an iterator that cannot be incremented
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T> >(grid_input,I)
    {}

    static RANGE<VECTOR<int,1> > Sentinels()
    {return RANGE<VECTOR<int,1> >(-1,0);}

private:
    void Initialize(const RANGE<VECTOR<int,1> >& region)
    {int index=grid.block_run_offsets(region.min_corner.x),index_end=grid.block_run_offsets(region.max_corner.x+1)-1;
    if(index<=index_end) RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T> >::Initialize(index,index_end);
    else{run=run_end=0;j_end=0;j=1;dj=0;}}
public:

    void Initialize(const TV& X)
    {grid.Clamped_Block(X,j,run);if(!run){j_end=0;j=1;dj=0;return;}
    j_end=j;dj=j-run->jmin;}

    void Initialize(const RLE_GRID_ITERATOR_BLOCK_2D<T>& block) // creates an iterator that cannot be incremented
    {RLE_GRID_ITERATOR_BLOCK<RLE_GRID_2D<T> >::Initialize(block);}

    void operator++(int)
    {dj++;
    if(++j>j_end){
        if(++run>run_end) return; // finished
        j_end=run->jmax;j=run->jmin;dj=0;}}

    VECTOR<int,2> I() const
    {return VECTOR<int,2>(run->i,j);}

    TV Minimum_Corner() const
    {return grid.uniform_grid.Center(run->i,j);}

    TV Center() const
    {return grid.uniform_grid.X(run->i+1,j+1);}

    RANGE<TV> Bounding_Box() const
    {TV minimum_corner=Minimum_Corner(),maximum_corner=minimum_corner+grid.uniform_grid.dX;
    return RANGE<TV>(minimum_corner,maximum_corner);}

    TV Cell_X(const int index) const
    {assert(*this && 0<=index&&index<4);TV X=Minimum_Corner();
    if(index&1) X.x+=grid.uniform_grid.dX.x;if(index&2) X.y+=grid.uniform_grid.dX.y;
    return X;}

    // TODO: these are not the typical parameters to Inside functions
    bool Inside(const TV& X,const T thickness_multiplier=-(T)1e-3) const
    {return Bounding_Box().Inside(X,thickness_multiplier*grid.Minimum_Edge_Length());}

    bool Lazy_Inside(const TV& X) const
    {return Bounding_Box().Lazy_Inside(X);}

//#####################################################################
};
}
#endif
#endif
