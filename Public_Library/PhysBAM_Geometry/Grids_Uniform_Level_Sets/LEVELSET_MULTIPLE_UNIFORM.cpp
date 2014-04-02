//#####################################################################
// Copyright 2005, Ron Fedkiw, Eran Guendelman, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MULTIPLE_UNIFORM  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_1D.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_2D.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//##################################################################### 
template<class T_GRID> LEVELSET_MULTIPLE_UNIFORM<T_GRID>::
LEVELSET_MULTIPLE_UNIFORM(T_GRID& grid_input,ARRAY<T_ARRAYS_SCALAR>& phis_input,const bool use_external_levelsets_input)
    :LEVELSET_MULTIPLE<T_GRID>(grid_input,phis_input,use_external_levelsets_input)
{}
//#####################################################################
// Destructor
//##################################################################### 
template<class T_GRID> LEVELSET_MULTIPLE_UNIFORM<T_GRID>::
~LEVELSET_MULTIPLE_UNIFORM()
{}
//#####################################################################
// Function Project_Levelset
//##################################################################### 
template<class T_GRID> void LEVELSET_MULTIPLE_UNIFORM<T_GRID>::
Project_Levelset(const int number_of_ghost_cells)
{
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
        int minimum_region,second_minimum_region;T minimum_phi,second_minimum_phi;
        Two_Minimum_Regions(iterator.Cell_Index(),minimum_region,second_minimum_region,minimum_phi,second_minimum_phi);
        T correction=(T).5*(minimum_phi+second_minimum_phi);
        for(int k=1;k<=phis.m;k++) phis(k)(iterator.Cell_Index())-=correction;}
}
//#####################################################################
// Function Project_Levelset
//##################################################################### 
template<class T_GRID> void LEVELSET_MULTIPLE_UNIFORM<T_GRID>::
Get_Single_Levelset(const ARRAY<bool>& positive_regions,T_LEVELSET& levelset,const bool flood_fill_for_bubbles)
{
    T_ARRAYS_SCALAR& phi_ghost=levelset.phi;phi_ghost.Resize(grid.Domain_Indices(3));
    if(flood_fill_for_bubbles){
        T_ARRAYS_INT colors(grid.Domain_Indices(3));T_FACE_ARRAYS_BOOL edge_is_blocked(grid,3);
        for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){
            if(!positive_regions(Inside_Region(iterator.Cell_Index()))) colors(iterator.Cell_Index())=-1;}
        T_FLOOD_FILL flood_fill;int number_of_colors=flood_fill.Flood_Fill(colors,edge_is_blocked);
        ARRAY<bool> color_touches_top_of_domain(number_of_colors);
        for(FACE_ITERATOR iterator(grid,0,T_GRID::BOUNDARY_REGION,4);iterator.Valid();iterator.Next()){
            if(colors(iterator.First_Cell_Index())>0)color_touches_top_of_domain(colors(iterator.First_Cell_Index()))=true;}
        for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){
            T minimum_phi;Inside_Region(iterator.Cell_Index(),minimum_phi);
            if(colors(iterator.Cell_Index())>0 && color_touches_top_of_domain(colors(iterator.Cell_Index())))
                phi_ghost(iterator.Cell_Index())=-minimum_phi; // make levelset positive in the dirichlet regions
            else phi_ghost(iterator.Cell_Index())=minimum_phi;}}
    else{
        for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){
            T minimum_phi;int minimum_region=Inside_Region(iterator.Cell_Index(),minimum_phi);
            if(positive_regions(minimum_region)) phi_ghost(iterator.Cell_Index())=-minimum_phi; // make levelset positive in the dirichlet regions
            else phi_ghost(iterator.Cell_Index())=minimum_phi;}}
    int number_of_positive_regions=0;for(int i=1;i<=positive_regions.m;i++)if(positive_regions(i))number_of_positive_regions++;
    if(number_of_positive_regions>1)levelset.Fast_Marching_Method();
}
//#####################################################################
template class LEVELSET_MULTIPLE_UNIFORM<GRID<VECTOR<float,1> > >;
template class LEVELSET_MULTIPLE_UNIFORM<GRID<VECTOR<float,2> > >;
template class LEVELSET_MULTIPLE_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_MULTIPLE_UNIFORM<GRID<VECTOR<double,1> > >;
template class LEVELSET_MULTIPLE_UNIFORM<GRID<VECTOR<double,2> > >;
template class LEVELSET_MULTIPLE_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
