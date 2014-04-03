//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR_FACE
//#####################################################################
#ifndef __UNIFORM_GRID_ITERATOR_FACE__
#define __UNIFORM_GRID_ITERATOR_FACE__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>

namespace PhysBAM{

template<class TV>
class UNIFORM_GRID_ITERATOR_FACE:public UNIFORM_GRID_ITERATOR<TV>
{
public:
    typedef typename GRID<TV>::REGION T_REGION;typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    typedef TV VECTOR_T;
    using UNIFORM_GRID_ITERATOR<TV>::grid;using UNIFORM_GRID_ITERATOR<TV>::index;using UNIFORM_GRID_ITERATOR<TV>::region;using UNIFORM_GRID_ITERATOR<TV>::valid;
    using UNIFORM_GRID_ITERATOR<TV>::Reset;using UNIFORM_GRID_ITERATOR<TV>::current_region;using UNIFORM_GRID_ITERATOR<TV>::Add_Region;
    using UNIFORM_GRID_ITERATOR<TV>::Reset_Regions;

protected:
    T_REGION region_type;
    int side;
    int axis;
    bool single_axis;
    int number_of_ghost_cells;
    T face_size;

public:
    // axis_input==0 means iterate through faces in all dimensions
    UNIFORM_GRID_ITERATOR_FACE(const GRID<TV>& grid_input,const int number_of_ghost_cells_input=0,const T_REGION& region_type_input=GRID<TV>::WHOLE_REGION,const int side_input=0,
        int axis_input=0);

    UNIFORM_GRID_ITERATOR_FACE(const GRID<TV>& grid_input,const int axis_input,const TV_INT& face_index);

    UNIFORM_GRID_ITERATOR_FACE(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,const int axis_input);

private:
    void Reset_Axis(const int axis_input);
    void Next_Helper();

public:
    void Next() PHYSBAM_ALWAYS_INLINE // overloads UNIFORM_GRID_ITERATOR::Next but we don't want that to be virtual to avoid virtual call overhead
    {if(index(TV::dimension)<region.max_corner(TV::dimension)) index(TV::dimension)++;else Next_Helper();}

    int Axis() const
    {return axis;}

    const TV_INT& Face_Index() const
    {return index;}

    FACE_INDEX<TV::dimension> Full_Index() const
    {return FACE_INDEX<TV::dimension>(axis,index);}

    T Face_Size() const
    {return face_size;}

    TV Location() const
    {return grid.Face(axis,index);}

    TV_INT First_Cell_Index() const
    {TV_INT i(index);i(axis)--;return i;}

    TV_INT Second_Cell_Index() const
    {return index;}

    void Unordered_Cell_Indices_Touching_Face(TV_INT& cell1,TV_INT& cell2) const
    {cell1=First_Cell_Index();cell2=Second_Cell_Index();}

    RANGE<TV> Dual_Cell() const
    {return RANGE<TV>(Location()-(T).5*grid.dX,Location()+(T).5*grid.dX);}

    TV First_Cell_Center() const
    {return grid.Center(First_Cell_Index());}

    TV Second_Cell_Center() const
    {return grid.Center(Second_Cell_Index());}

    bool First_Boundary() const // returns true if currently on left, bottom, or front boundary
    {assert(region_type==GRID<TV>::BOUNDARY_REGION);return (!side && current_region%2==0) || (side && (side-1)%2==0);}

    TV_INT Face_Node_Index(const int node) const // 1-based
    {return grid.Face_Node_Index(axis,index,node);}
//#####################################################################
};
}
#endif
