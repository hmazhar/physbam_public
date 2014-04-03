//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR_EDGE
//#####################################################################
#ifndef __UNIFORM_GRID_ITERATOR_EDGE__
#define __UNIFORM_GRID_ITERATOR_EDGE__

#include <PhysBAM_Tools/Grids_Uniform/EDGE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>

namespace PhysBAM{

template<class TV>
class UNIFORM_GRID_ITERATOR_EDGE:public UNIFORM_GRID_ITERATOR<TV>
{
public:
    typedef typename GRID<TV>::REGION T_REGION;typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    using UNIFORM_GRID_ITERATOR<TV>::grid;using UNIFORM_GRID_ITERATOR<TV>::index;using UNIFORM_GRID_ITERATOR<TV>::region;using UNIFORM_GRID_ITERATOR<TV>::valid;
    using UNIFORM_GRID_ITERATOR<TV>::Reset;using UNIFORM_GRID_ITERATOR<TV>::current_region;using UNIFORM_GRID_ITERATOR<TV>::Add_Region;
    using UNIFORM_GRID_ITERATOR<TV>::Reset_Regions;

protected:
    T_REGION region_type;
    int side;
    int axis;
    bool single_axis;
    int number_of_ghost_cells;

public:
    // axis_input==0 means iterate through faces in all dimensions
    UNIFORM_GRID_ITERATOR_EDGE(const GRID<TV>& grid_input,const int number_of_ghost_cells_input=0,const T_REGION& region_type_input=GRID<TV>::WHOLE_REGION,const int side_input=0,
        int axis_input=0);

    UNIFORM_GRID_ITERATOR_EDGE(const GRID<TV>& grid_input,const int axis_input,const TV_INT& face_index);

    UNIFORM_GRID_ITERATOR_EDGE(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,const int axis_input);

private:
    void Reset_Axis(const int axis_input);
    void Next_Helper();

public:
    void Next() PHYSBAM_ALWAYS_INLINE // overloads UNIFORM_GRID_ITERATOR::Next but we don't want that to be virtual to avoid virtual call overhead
    {if(index(TV::dimension)<region.max_corner(TV::dimension)) index(TV::dimension)++;else Next_Helper();}

    int Axis() const
    {return axis;}

    EDGE_INDEX<TV::dimension> Full_Index() const
    {return EDGE_INDEX<TV::dimension>(axis,index);}
//#####################################################################
};
}
#endif
