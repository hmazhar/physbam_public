//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR_EDGE
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_EDGE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
// axis_input==0 means iterate through faces in all dimensions
template<class TV> UNIFORM_GRID_ITERATOR_EDGE<TV>::
UNIFORM_GRID_ITERATOR_EDGE(const GRID<TV>& grid_input,const int number_of_ghost_cells_input,const T_REGION& region_type_input,const int side_input,int axis_input)
    :UNIFORM_GRID_ITERATOR<TV>(grid_input),region_type(region_type_input),side(side_input),number_of_ghost_cells(number_of_ghost_cells_input)
{
    assert(side==0);
    assert(region_type==GRID<TV>::WHOLE_REGION);
    if(axis_input){single_axis=true;Reset_Axis(axis_input);}
    else{single_axis=false;Reset_Axis(1);}
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_GRID_ITERATOR_EDGE<TV>::
UNIFORM_GRID_ITERATOR_EDGE(const GRID<TV>& grid_input,const int axis_input,const TV_INT& face_index)
    :UNIFORM_GRID_ITERATOR<TV>(grid_input,RANGE<TV_INT>(face_index))
{
    single_axis=true;
    axis=axis_input;
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_GRID_ITERATOR_EDGE<TV>::
UNIFORM_GRID_ITERATOR_EDGE(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,const int axis_input)
    :UNIFORM_GRID_ITERATOR<TV>(grid_input),axis(axis_input),single_axis(true)
{
    assert(axis);
    Add_Region(explicit_region_input);
    Reset();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void UNIFORM_GRID_ITERATOR_EDGE<TV>::
Next_Helper()
{
    UNIFORM_GRID_ITERATOR<TV>::Next_Helper();
    if(!valid && !single_axis && axis<TV::dimension) Reset_Axis(axis+1);
}
//#####################################################################
// Function Reset_Axis
//#####################################################################
template<class TV> void UNIFORM_GRID_ITERATOR_EDGE<TV>::
Reset_Axis(const int axis_input)
{
    axis=axis_input;Reset_Regions();
    RANGE<TV_INT> domain(grid.Node_Indices(number_of_ghost_cells));
    switch(region_type){
        case GRID<TV>::WHOLE_REGION:
            assert(!side);
            domain.max_corner(axis)--;
            Add_Region(domain);
            break;
        default:PHYSBAM_NOT_IMPLEMENTED();}
    Reset();
}
//#####################################################################
template class UNIFORM_GRID_ITERATOR_EDGE<VECTOR<float,1> >;
template class UNIFORM_GRID_ITERATOR_EDGE<VECTOR<float,2> >;
template class UNIFORM_GRID_ITERATOR_EDGE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_GRID_ITERATOR_EDGE<VECTOR<double,1> >;
template class UNIFORM_GRID_ITERATOR_EDGE<VECTOR<double,2> >;
template class UNIFORM_GRID_ITERATOR_EDGE<VECTOR<double,3> >;
#endif
