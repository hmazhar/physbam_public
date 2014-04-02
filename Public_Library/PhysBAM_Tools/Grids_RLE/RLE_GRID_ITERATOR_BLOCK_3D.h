//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_BLOCK_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATOR_BLOCK_3D__
#define __RLE_GRID_ITERATOR_BLOCK_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK.h>
namespace PhysBAM{

template<class T_input>
class RLE_GRID_ITERATOR_BLOCK_3D:public RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T_input> >
{
    typedef T_input T;
public:
    typedef RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T> > BASE;
    using BASE::grid;using BASE::run_end;using BASE::run;using BASE::j_end;using BASE::j;using BASE::dj;using BASE::Face_X;using BASE::Face_Y;using BASE::Face_Z;

    int i,i_end,ij_start,ij_end;

    RLE_GRID_ITERATOR_BLOCK_3D(const RLE_GRID_3D<T>& grid_input,const int number_of_ghost_cells)
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T> >(grid_input)
    {
        RANGE<VECTOR<int,2> > region=grid.horizontal_grid.Get_MAC_Grid().Domain_Indices()+Sentinels();
        Initialize(region.Thickened(number_of_ghost_cells-1));
    }

    RLE_GRID_ITERATOR_BLOCK_3D(const RLE_GRID_3D<T>& grid_input,const RANGE<VECTOR<int,2> >& region)
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T> >(grid_input)
    {
        Initialize(region);
    }

    RLE_GRID_ITERATOR_BLOCK_3D(const RLE_GRID_3D<T>& grid_input,const VECTOR<T,3>& X) // creates an iterator that cannot be incremented
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T> >(grid_input)
    {
        Initialize(X);
    }

    RLE_GRID_ITERATOR_BLOCK_3D(const RLE_GRID_3D<T>& grid_input,const VECTOR<int,3>& I) // creates an iterator that cannot be incremented
        :RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T> >(grid_input,I),i(I.x)
    {}

    static RANGE<VECTOR<int,2> > Sentinels()
    {return RANGE<VECTOR<int,2> >(-1,0,-1,0);}

private:
    void Initialize(const RANGE<VECTOR<int,2> >& region)
    {i=region.min_corner.x;i_end=region.max_corner.x;ij_start=region.min_corner.y;ij_end=region.max_corner.y;Initialize_Slice();}

    void Initialize_Slice()
    {for(;i<=i_end;i++){
        int index=grid.block_run_offsets(i,ij_start),index_end=grid.block_run_offsets(i,ij_end+1)-1;
        if(index<=index_end){RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T> >::Initialize(index,index_end);return;}}
    run=run_end=0;j_end=0;j=1;dj=0;} // finished
public:

    void Initialize(const VECTOR<T,3>& X)
    {grid.Clamped_Block(X,i,j,run);if(!run){j_end=0;j=1;dj=0;return;}
    j_end=j;dj=j-run->jmin;}

    void Initialize(const RLE_GRID_ITERATOR_BLOCK_3D<T>& block) // creates an iterator that cannot be incremented
    {RLE_GRID_ITERATOR_BLOCK<RLE_GRID_3D<T> >::Initialize(block);i=block.i;}

    void operator++(int)
    {dj++;
    if(++j>j_end){
        if(++run>run_end){i++;Initialize_Slice();}
        else{j_end=run->jmax;j=run->jmin;dj=0;}}}

    VECTOR<int,3> I() const
    {return VECTOR<int,3>(i,j,run->ij);}

    VECTOR<T,3> Minimum_Corner() const
    {return grid.uniform_grid.Center(i,j,run->ij);}

    RANGE<VECTOR<T,3> > Maximum_Corner() const
    {return grid.uniform_grid.Center(run->i+1,j+1,run->ij+1);}

    VECTOR<T,3> Center() const
    {return grid.uniform_grid.X(i+1,j+1,run->ij+1);}

    RANGE<VECTOR<T,3> > Bounding_Box() const
    {VECTOR<T,3> minimum_corner=Minimum_Corner(),maximum_corner=minimum_corner+grid.uniform_grid.dX;
    return RANGE<VECTOR<T,3> >(minimum_corner,maximum_corner);}

    VECTOR<T,3> Cell_X(const int index) const
    {assert(*this && 0<=index&&index<8);VECTOR<T,3> X=Minimum_Corner();
    if(index&1) X.x+=grid.uniform_grid.dX.x;if(index&2) X.y+=grid.uniform_grid.dX.y;if(index&4) X.z+=grid.uniform_grid.dX.z;
    return X;}

    // TODO: these are not the typical parameters to Inside functions
    bool Inside(const VECTOR<T,3>& X,const T thickness_multiplier=-(T)1e-3) const
    {return Bounding_Box().Inside(X,thickness_multiplier*grid.Minimum_Edge_Length());}

    bool Lazy_Inside(const VECTOR<T,3>& X) const
    {return Bounding_Box().Lazy_Inside(X);}

//#####################################################################
};
}
#endif
#endif
