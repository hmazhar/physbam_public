//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace  GRADIENT_UNIFORM
//##################################################################### 
#ifndef __GRADIENT_UNIFORM__
#define __GRADIENT_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

namespace GRADIENT{

template<class T,int d>
void Compute_Magnitude(const GRID<VECTOR<T,d> >& grid,const int number_of_ghost_cells,ARRAY<T,VECTOR<int,d> >& values,ARRAY<T,VECTOR<int,d> >& gradient){
    typedef typename GRID<VECTOR<T,d> >::CELL_ITERATOR CELL_ITERATOR;

    VECTOR<T,d> one_over_dx=grid.one_over_dX;
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells-1);iterator.Valid();iterator.Next()){
        VECTOR<int,d> cell_index=iterator.Cell_Index();T sum_of_partials=0;
        for(int axis=1;axis<=d;axis++){
            T partial=(values(iterator.Cell_Neighbor(2*axis))-values(iterator.Cell_Neighbor(2*axis-1)))*(T).5*one_over_dx[axis];
            sum_of_partials+=partial*partial;}
        gradient(cell_index)=sqrt(sum_of_partials);}
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells,GRID<VECTOR<T,d> >::BOUNDARY_INTERIOR_REGION);iterator.Valid();iterator.Next()) gradient(iterator.Cell_Index())=(T)0;
}
}
}
#endif
