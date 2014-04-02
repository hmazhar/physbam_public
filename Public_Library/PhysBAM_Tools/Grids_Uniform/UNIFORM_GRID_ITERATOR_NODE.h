//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR_NODE
//#####################################################################
#ifndef __UNIFORM_GRID_ITERATOR_NODE__
#define __UNIFORM_GRID_ITERATOR_NODE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>
namespace PhysBAM{

template<class TV>
class UNIFORM_GRID_ITERATOR_NODE:public UNIFORM_GRID_ITERATOR<TV>
{
public:
    typedef typename GRID<TV>::REGION T_REGION;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef TV VECTOR_T;
    using UNIFORM_GRID_ITERATOR<TV>::grid;using UNIFORM_GRID_ITERATOR<TV>::index;using UNIFORM_GRID_ITERATOR<TV>::Add_Region;using UNIFORM_GRID_ITERATOR<TV>::Reset;

    UNIFORM_GRID_ITERATOR_NODE(const GRID<TV>& grid_input,const int number_of_ghost_cells=0,const T_REGION& region_type=GRID<TV>::WHOLE_REGION,const int side=0);
    UNIFORM_GRID_ITERATOR_NODE(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input);

    const TV_INT& Node_Index() const
    {return index;}

    TV Location() const
    {return grid.Node(index);}

    TV_INT Node_Neighbor(const int i) const // i=1 to 6
    {return grid.Node_Neighbor(index,i);}

//#####################################################################
};
}
#endif
