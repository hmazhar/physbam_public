//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_GRID_ITERATOR_CELL
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
using namespace PhysBAM;
template<class TV> UNIFORM_GRID_ITERATOR_CELL<TV>::
UNIFORM_GRID_ITERATOR_CELL(const GRID<TV>& grid_input,const int number_of_ghost_cells,const T_REGION& region_type,const int side)
    :UNIFORM_GRID_ITERATOR<TV>(grid_input)
{
    assert(0<=side&&side<=6);assert(region_type!=GRID<TV>::BOUNDARY_REGION && region_type!=GRID<TV>::INTERIOR_REGION); // these types may not really make sense here
    RANGE<TV_INT> domain(grid.Domain_Indices(number_of_ghost_cells));
    switch(region_type){
        case GRID<TV>::WHOLE_REGION: assert(!side);Add_Region(domain);break;
        case GRID<TV>::BOUNDARY_INTERIOR_REGION: // outer boundary of grid with specified ghost cells
            if(!side){ // don't loop over the same node twice!
                RANGE<TV_INT> domain_copy(domain);
                for(int axis=TV::dimension;axis>=1;axis--){
                    domain.max_corner(axis)=domain.min_corner(axis);
                    Add_Region(domain);
                    domain.max_corner(axis)=domain.min_corner(axis)=domain_copy.max_corner(axis);
                    Add_Region(domain);
                    domain.max_corner(axis)=domain_copy.max_corner(axis)-1;
                    domain.min_corner(axis)=domain_copy.min_corner(axis)+1;}}
            else{int axis=(side+1)/2;if(side&1) domain.max_corner(axis)=domain.min_corner(axis);else domain.min_corner(axis)=domain.max_corner(axis);Add_Region(domain);}
            break;
        default: assert(region_type==GRID<TV>::GHOST_REGION && number_of_ghost_cells>0); // ghost region of grid with specified ghost cells
            if(!side){ // don't loop over the same cell twice!
                TV_INT max_copy(domain.max_corner);
                for(int axis=TV::dimension;axis>=1;axis--){
                    domain.max_corner(axis)=0;
                    Add_Region(domain);
                    domain.max_corner(axis)=max_copy(axis);
                    domain.min_corner(axis)=grid.numbers_of_cells(axis)+1;
                    Add_Region(domain);
                    domain.max_corner(axis)=grid.numbers_of_cells(axis);
                    domain.min_corner(axis)=1;}}
            else{int axis=(side+1)/2;if(side&1) domain.max_corner(axis)=0;else domain.min_corner(axis)=grid.numbers_of_cells(axis)+1;Add_Region(domain);}
            break;}
    Reset();
}
//#####################################################################
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<float,0> >;
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<float,1> >;
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<float,2> >;
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<double,0> >;
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<double,1> >;
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<double,2> >;
template class UNIFORM_GRID_ITERATOR_CELL<VECTOR<double,3> >;
#endif
