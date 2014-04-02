//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR_FACE
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
// axis_input==0 means iterate through faces in all dimensions
template<class TV> UNIFORM_GRID_ITERATOR_FACE<TV>::
UNIFORM_GRID_ITERATOR_FACE(const GRID<TV>& grid_input,const int number_of_ghost_cells_input,const T_REGION& region_type_input,const int side_input,int axis_input)
    :UNIFORM_GRID_ITERATOR<TV>(grid_input),region_type(region_type_input),side(side_input),number_of_ghost_cells(number_of_ghost_cells_input)
{
    assert(0<=side&&side<=2*TV::dimension&&(!side||!axis_input||(side+1)/2==axis_input));
    assert(region_type!=GRID<TV>::BOUNDARY_INTERIOR_REGION); // TODO: implement this case!
    if(region_type==GRID<TV>::BOUNDARY_REGION && side) axis_input=(side+1)/2;
    if(axis_input){single_axis=true;Reset_Axis(axis_input);}else{single_axis=false;Reset_Axis(1);}
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_GRID_ITERATOR_FACE<TV>::
UNIFORM_GRID_ITERATOR_FACE(const GRID<TV>& grid_input,const int axis_input,const TV_INT& face_index)
    :UNIFORM_GRID_ITERATOR<TV>(grid_input,RANGE<TV_INT>(face_index))
{
    single_axis=true;axis=axis_input;face_size=grid.Face_Size(axis);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_GRID_ITERATOR_FACE<TV>::
UNIFORM_GRID_ITERATOR_FACE(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,const int axis_input)
    :UNIFORM_GRID_ITERATOR<TV>(grid_input),axis(axis_input),single_axis(true)
{
    assert(axis);Add_Region(explicit_region_input);face_size=grid.Face_Size(axis);Reset();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void UNIFORM_GRID_ITERATOR_FACE<TV>::
Next_Helper()
{
    UNIFORM_GRID_ITERATOR<TV>::Next_Helper();
    if(!valid && !single_axis && axis<TV::dimension) Reset_Axis(axis+1);
}
//#####################################################################
// Function Reset_Axis
//#####################################################################
template<class TV> void UNIFORM_GRID_ITERATOR_FACE<TV>::
Reset_Axis(const int axis_input)
{
    axis=axis_input;Reset_Regions();
    RANGE<TV_INT> domain(grid.Cell_Indices(number_of_ghost_cells));
    switch(region_type){
        case GRID<TV>::WHOLE_REGION:{assert(!side);
            domain.max_corner(axis)++;
            Add_Region(domain);
            break;}
        case GRID<TV>::GHOST_REGION:{
            if(!side){ // TODO(jontg): Beware of duplicates!
                for(int side_iterator=1;side_iterator<=TV::dimension*2;++side_iterator){
                    int axis_of_side=(side_iterator+1)/2;
                    if(side_iterator&1){
                        RANGE<TV_INT> domain_copy(domain); domain_copy.max_corner(axis)++;
                        domain_copy.max_corner(axis_of_side)=domain_copy.min_corner(axis_of_side)+number_of_ghost_cells-1;Add_Region(domain_copy);}
                    if(!(side_iterator&1)){
                        RANGE<TV_INT> domain_copy(domain); domain_copy.max_corner(axis)++;
                        domain_copy.min_corner(axis_of_side)=domain_copy.max_corner(axis_of_side)-number_of_ghost_cells+1;Add_Region(domain_copy);}}}
            else{
                int axis_of_side=(side+1)/2;
                if(side&1){
                    RANGE<TV_INT> domain_copy(domain); domain_copy.max_corner(axis)++;
                    domain_copy.max_corner(axis_of_side)=domain_copy.min_corner(axis_of_side)+number_of_ghost_cells-1;Add_Region(domain_copy);}
                if(!(side&1)){
                    RANGE<TV_INT> domain_copy(domain); domain_copy.max_corner(axis)++;
                    domain_copy.min_corner(axis_of_side)=domain_copy.max_corner(axis_of_side)-number_of_ghost_cells+1;Add_Region(domain_copy);}}
            break;}
        case GRID<TV>::BOUNDARY_REGION:{
            if(!side || side&1){RANGE<TV_INT> domain_copy(domain);domain_copy.max_corner(axis)=domain.min_corner(axis);Add_Region(domain_copy);}
            if(!side || !(side&1)){RANGE<TV_INT> domain_copy(domain);domain_copy.min_corner(axis)=domain_copy.max_corner(axis)=domain.max_corner(axis)+1;Add_Region(domain_copy);}
            break;}
        default:{assert(region_type==GRID<TV>::INTERIOR_REGION && !side);
            domain.min_corner(axis)++;
            Add_Region(domain);
            break;}}
    face_size=grid.Face_Size(axis);Reset();
}
//#####################################################################
template class UNIFORM_GRID_ITERATOR_FACE<VECTOR<float,1> >;
template class UNIFORM_GRID_ITERATOR_FACE<VECTOR<float,2> >;
template class UNIFORM_GRID_ITERATOR_FACE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_GRID_ITERATOR_FACE<VECTOR<double,1> >;
template class UNIFORM_GRID_ITERATOR_FACE<VECTOR<double,2> >;
template class UNIFORM_GRID_ITERATOR_FACE<VECTOR<double,3> >;
#endif
