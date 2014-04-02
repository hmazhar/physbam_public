//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR_CELL
//#####################################################################
#ifndef __UNIFORM_GRID_ITERATOR_CELL__
#define __UNIFORM_GRID_ITERATOR_CELL__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>
namespace PhysBAM{

template<class TV>
class UNIFORM_GRID_ITERATOR_CELL:public UNIFORM_GRID_ITERATOR<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;
public:
    typedef typename GRID<TV>::REGION T_REGION;
    typedef TV VECTOR_T;
    using UNIFORM_GRID_ITERATOR<TV>::grid;using UNIFORM_GRID_ITERATOR<TV>::index;using UNIFORM_GRID_ITERATOR<TV>::Add_Region;using UNIFORM_GRID_ITERATOR<TV>::Reset;

    UNIFORM_GRID_ITERATOR_CELL(const GRID<TV>& grid_input,const int number_of_ghost_cells=0,const T_REGION& region_type=GRID<TV>::WHOLE_REGION,const int side=0);

    UNIFORM_GRID_ITERATOR_CELL(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input)
        :UNIFORM_GRID_ITERATOR<TV>(grid_input,region_input)
    {}

    const TV_INT& Cell_Index() const
    {return index;}

    TV Location() const
    {return grid.Center(index);}

    RANGE<TV> Bounding_Box() const
    {TV minimum_corner=grid.Node(index);return RANGE<TV>(minimum_corner,minimum_corner+grid.dX);}

    TV_INT Cell_Node_Index(const int node) const
    {return index+GRID<TV>::Binary_Counts(index)(node);}

    TV_INT First_Face_Index(const int axis) const
    {return grid.First_Face_Index_In_Cell(axis,index);}

    TV_INT Second_Face_Index(const int axis) const
    {return grid.Second_Face_Index_In_Cell(axis,index);}

    FACE_INDEX<TV::dimension> Full_First_Face_Index(const int axis) const
    {return FACE_INDEX<TV::dimension>(axis,grid.First_Face_Index_In_Cell(axis,index));}

    FACE_INDEX<TV::dimension> Full_Second_Face_Index(const int axis) const
    {return FACE_INDEX<TV::dimension>(axis,grid.Second_Face_Index_In_Cell(axis,index));}

    TV_INT Cell_Neighbor(const int i) const // 1 to 2*dimension
    {return GRID<TV>::Node_Neighbor(index,i);}
};
//#####################################################################
}
#endif
