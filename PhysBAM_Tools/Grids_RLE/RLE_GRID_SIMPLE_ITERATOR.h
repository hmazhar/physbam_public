//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_SIMPLE_ITERATOR
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_SIMPLE_ITERATOR__
#define __RLE_GRID_SIMPLE_ITERATOR__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
namespace PhysBAM{

//#####################################################################
// Class RLE_GRID_SIMPLE_ITERATOR<T_GRID,T_ITERATOR>
//#####################################################################
template<class T_GRID,class T_ITERATOR>
class RLE_GRID_SIMPLE_ITERATOR
{
    RLE_GRID_SIMPLE_ITERATOR(); // private constructor, since only dimension specific template specializations should be used
};
//#####################################################################
// Class RLE_GRID_SIMPLE_ITERATOR<RLE_GRID_2D<T>,T_ITERATOR>
//#####################################################################
template<class T,class T_ITERATOR>
class RLE_GRID_SIMPLE_ITERATOR<RLE_GRID_2D<T>,T_ITERATOR>
{
    int index_end;
public:
    int index;

    RLE_GRID_SIMPLE_ITERATOR(const RLE_GRID_2D<T>& grid,const RANGE<VECTOR<int,1> >& region)
    {
        RANGE<VECTOR<int,1> > indices=T_ITERATOR::Indices_In_Slices(grid,region.min_corner.x,region.max_corner.x);
        index=indices.min_corner.x;index_end=indices.max_corner.x;
    }

    operator bool() const
    {return index<=index_end;}

    void operator++(int)
    {assert(*this);index++;}
};
//#####################################################################
// Class RLE_GRID_SIMPLE_ITERATOR<RLE_GRID_3D<T>,T_ITERATOR>
//#####################################################################
template<class T,class T_ITERATOR>
class RLE_GRID_SIMPLE_ITERATOR<RLE_GRID_3D<T>,T_ITERATOR>
{
    const RLE_GRID_3D<T>& grid;
    const int i_end,ij_start,ij_end;
    int i,index_end;
public:
    int index;

    RLE_GRID_SIMPLE_ITERATOR(const RLE_GRID_3D<T>& grid_input,const RANGE<VECTOR<int,2> >& region)
        :grid(grid_input),i_end(region.max_corner.x),ij_start(region.min_corner.y),ij_end(region.max_corner.y),i(region.min_corner.x)
    {
        Start_Slice();
    }

private:
    void Start_Slice()
    {RANGE<VECTOR<int,1> > indices=T_ITERATOR::Indices_In_Slice(grid,i,ij_start,ij_end);index=indices.min_corner.x;index_end=indices.max_corner.x;}
public:

    operator bool() const
    {return i<=i_end;}

    void operator++(int)
    {assert(*this);
    if(++index>index_end){
        if(++i>i_end) return; // finished
        Start_Slice();}}
};
//#####################################################################
}
#endif
#endif
