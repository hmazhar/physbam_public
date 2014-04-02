//#####################################################################
// Copyright 2010, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace  SMOOTH_UNIFORM
//##################################################################### 
#ifndef __SMOOTH_UNIFORM__
#define __SMOOTH_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>

namespace PhysBAM{

namespace SMOOTH{

//#####################################################################
// Function Smooth
//#####################################################################
// explicit update with CFL=1/2 and all Neumann outer boundaries, no smoothing where phi is negative if defined
template<class T_GRID,class T_ARRAYS_T2>
void Smooth(T_ARRAYS_T2& d,const int steps,const typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR* phi){
    typedef typename T_ARRAYS_T2::ELEMENT T2;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::SCALAR T;

    if(!steps) return;
    int number_of_ghost_cells=1-d.domain.min_corner.x;
    RANGE<TV_INT> domain_indices=d.Domain_Indices();
    T_GRID grid(domain_indices.Thickened(-number_of_ghost_cells).Maximum_Corner(),RANGE<typename T_GRID::VECTOR_T>::Centered_Box());
    if(grid.Domain_Indices().Thickened(number_of_ghost_cells)!=d.Domain_Indices()) PHYSBAM_FATAL_ERROR();

    T_ARRAYS_T2 d_ghost(grid.Domain_Indices(number_of_ghost_cells+1),false);    
    for(int step=1;step<=steps;step++){
        T_ARRAYS_T2::Put(d,d_ghost);
        for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){
            int side=2*(axis-1)+axis_side;TV_INT ghost_offset=(axis_side==1?-1:1)*TV_INT::Axis_Vector(axis);
            for(NODE_ITERATOR iterator(grid,number_of_ghost_cells,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
                d_ghost(node+ghost_offset)=d_ghost(node);}}
        if(!phi) for(NODE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            T2 sum=(T)T_GRID::number_of_neighbors_per_node*d_ghost(node);
            for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++)sum+=d_ghost(T_GRID::Node_Neighbor(node,n));
            d(node)=(T)1/(2*T_GRID::number_of_neighbors_per_node)*sum;}
        else for(NODE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();if((*phi)(node) > 0){
            T2 sum=(T)T_GRID::number_of_neighbors_per_node*d_ghost(node);
            for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++)sum+=d_ghost(T_GRID::Node_Neighbor(node,n));
            d(node)=(T)1/(2*T_GRID::number_of_neighbors_per_node)*sum;}}}
}
}
}
#endif
