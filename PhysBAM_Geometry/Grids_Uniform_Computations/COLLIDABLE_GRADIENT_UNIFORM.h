//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace  COLLIDABLE_GRADIENT_UNIFORM
//##################################################################### 
#ifndef __COLLIDABLE_GRADIENT_UNIFORM__
#define __COLLIDABLE_GRADIENT_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
namespace PhysBAM{

namespace ARRAYS_COMPUTATIONS{

template<class T,int d>
void Collidable_Gradient_Magnitude(const typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,d> > >::GRID_BASED_COLLISION_GEOMETRY& collision_body_list,const GRID<VECTOR<T,d> >& grid,
                                   const int number_of_ghost_cells, ARRAY<T,VECTOR<int,d> >& values,ARRAY<T,VECTOR<int,d> >& gradient){
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
    typedef GRID<TV> T_GRID;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;

    T_FACE_ARRAYS_BOOL occluded_faces(grid,number_of_ghost_cells);
    collision_body_list.Compute_Psi_N(occluded_faces);
    for(FACE_ITERATOR iterator(grid,number_of_ghost_cells,T_GRID::BOUNDARY_INTERIOR_REGION);iterator.Valid();iterator.Next()) occluded_faces(iterator.Full_Index())=true;

    TV one_over_dx=grid.one_over_dX;
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells-1);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();T sum_of_partials=0;int stencil_width=0;
        for(int axis=1;axis<=d;axis++){T partial=0;
            bool left_cell_visible=!occluded_faces(iterator.Full_First_Face_Index(axis)),right_cell_visible=!occluded_faces(iterator.Full_Second_Face_Index(axis));
            if(left_cell_visible && right_cell_visible)
                partial=(values(iterator.Cell_Neighbor(2*axis))-values(iterator.Cell_Neighbor(2*axis-1)))*(T).5*one_over_dx[axis];
            else if(left_cell_visible)
                partial=(values(iterator.Cell_Index())-values(iterator.Cell_Neighbor(2*axis-1)))*one_over_dx[axis];
            else if(right_cell_visible)
                partial=(values(iterator.Cell_Neighbor(2*axis))-values(iterator.Cell_Index()))*one_over_dx[axis];
            else continue;
            sum_of_partials+=partial*partial;}
        gradient(cell_index)=sqrt(sum_of_partials);}
}
}
}
#endif
