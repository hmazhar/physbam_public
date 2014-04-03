//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING_METHOD_RLE  
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __FAST_MARCHING_METHOD_RLE__
#define __FAST_MARCHING_METHOD_RLE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
namespace PhysBAM{

template<class T_GRID>
class FAST_MARCHING_METHOD_RLE:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::BLOCK T_BLOCK;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
private:
    const LEVELSET_RLE<T_GRID>& levelset;
    const T_GRID& grid;
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& neighbors;
    const ARRAY<VECTOR<bool,T_GRID::dimension> >* neighbors_visible;
public:

    FAST_MARCHING_METHOD_RLE(const LEVELSET_RLE<T_GRID>& levelset_input);
    ~FAST_MARCHING_METHOD_RLE();

private:
    int Neighbor_If_Visible(const int k,const int cell) const
    {int neighbor=neighbors(cell)(k);return neighbor && (!neighbors_visible||(*neighbors_visible)(k&1?neighbor:cell)(((k-1)>>1)+1))?neighbor:0;}

    bool Neighbor_Visible(const int axis,const int cell1,const int cell2) const
    {assert(neighbors_visible);return (*neighbors_visible)(min(cell1,cell2))(axis);}
public:

//#####################################################################
    void Fast_Marching(ARRAY<T>& phi,const T stopping_distance,const bool outside_only=false) const;
    void Slow_Marching(ARRAY<T>& phi,const T stopping_distance,const T tolerance=(T).01) const;
private:
    void Update_Or_Add_Neighbor(ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,const int neighbor) const;
    void Initialize_Interface(ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,const bool outside_only=false) const;
    struct Find_Interface_Cells{template<class T_FACE> static void Apply(const FAST_MARCHING_METHOD_RLE<T_GRID>& fast_marching,ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done);};
    void Update_Close_Point(ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int cell) const;
//#####################################################################
};
}
#endif
#endif
