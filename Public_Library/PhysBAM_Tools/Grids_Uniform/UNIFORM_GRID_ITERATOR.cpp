//#####################################################################
// Copyright 2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_GRID_ITERATOR<TV>::
UNIFORM_GRID_ITERATOR(const GRID<TV>& grid_input)
    :grid(grid_input),number_of_regions(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_GRID_ITERATOR<TV>::
UNIFORM_GRID_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input)
    :grid(grid_input),number_of_regions(0)
{
    Add_Region(region_input);Reset();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void UNIFORM_GRID_ITERATOR<TV>::
Next_Helper()
{
    index(TV::dimension)=region.min_corner(TV::dimension);
    for(int i=TV::dimension-1;i>=1;i--){
        if(index(i)<region.max_corner(i)){index(i)++;return;}
        index(i)=region.min_corner(i);}
    Reset(current_region+1);
}
//#####################################################################
template class UNIFORM_GRID_ITERATOR<VECTOR<float,0> >;
template class UNIFORM_GRID_ITERATOR<VECTOR<float,1> >;
template class UNIFORM_GRID_ITERATOR<VECTOR<float,2> >;
template class UNIFORM_GRID_ITERATOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_GRID_ITERATOR<VECTOR<double,0> >;
template class UNIFORM_GRID_ITERATOR<VECTOR<double,1> >;
template class UNIFORM_GRID_ITERATOR<VECTOR<double,2> >;
template class UNIFORM_GRID_ITERATOR<VECTOR<double,3> >;
#endif
