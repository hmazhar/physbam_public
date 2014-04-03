//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Computations/VORTICITY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Function Vorticity
//#####################################################################
template<class T> void VORTICITY_DYADIC<T>::
Vorticity(QUADTREE_GRID<T>& grid,const ARRAY<VECTOR<T,2> >& V_ghost,ARRAY<T>& vorticity)
{
    T minimum_cell_size=grid.Minimum_Edge_Length();T one_over_two_dx=1/(2*minimum_cell_size);VECTOR<T,2> x_offset(minimum_cell_size,0),y_offset(0,minimum_cell_size);
    for(int i=1;i<=grid.number_of_nodes;i++)
        vorticity(i)=-(LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,grid.Node_Location(i)+y_offset).x
            -LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,grid.Node_Location(i)-y_offset).x)*one_over_two_dx;
    for(int i=1;i<=grid.number_of_nodes;i++)
        vorticity(i)+=(LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,grid.Node_Location(i)+x_offset).y
            -LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,grid.Node_Location(i)-x_offset).y)*one_over_two_dx;
}
//#####################################################################
// Function Vorticity
//#####################################################################
template<class T> void VORTICITY_DYADIC<T>::
Vorticity(OCTREE_GRID<T>& grid,const ARRAY<VECTOR<T,3> >& V_ghost,ARRAY<VECTOR<T,3> >& vorticity)
{
    T minimum_cell_size=grid.Minimum_Edge_Length();T one_over_two_dx=1/(2*minimum_cell_size);VECTOR<T,3> x_offset(minimum_cell_size,0,0),y_offset(0,minimum_cell_size,0),
        z_offset(0,0,minimum_cell_size);
    PHYSBAM_NOT_IMPLEMENTED(); // The following loop won't work since it touches the extreme boundary nodes, and the looks further out.  Also, it should probably be looping over cells.
    for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){
        int node=iterator.Node_Index();VECTOR<T,3> X=iterator.Location();
        vorticity(node).y=(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+z_offset).x
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-z_offset).x)*one_over_two_dx;
        vorticity(node).z=-(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+y_offset).x
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-y_offset).x)*one_over_two_dx;
        vorticity(node).x=-(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+z_offset).y
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-z_offset).y)*one_over_two_dx;
        vorticity(node).z+=(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+x_offset).y
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-x_offset).y)*one_over_two_dx;
        vorticity(node).x+=(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+y_offset).z
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-y_offset).z)*one_over_two_dx;
        vorticity(node).y+=-(LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X+x_offset).z
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,V_ghost,X-x_offset).z)*one_over_two_dx;}
}
//#####################################################################
template class VORTICITY_DYADIC<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTICITY_DYADIC<double>;
#endif
#endif
