//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_UNIFORM
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_REFLECTION_UNIFORM<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    T_ARRAYS_T2::Put(u,u_ghost); // interior
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++) Fill_Single_Ghost_Region(grid,u_ghost,regions(side),side,dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_REFLECTION_UNIFORM<T_GRID,T2>::
Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const
{
    if(Constant_Extrapolation(side)) Fill_Single_Ghost_Region(grid,u_ghost,side,region);
    else{int axis=(side+1)/2,boundary=Boundary(side,region),reflection_times_two=2*boundary+(side&1?-1:1);
        for(NODE_ITERATOR iterator(grid,region);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            TV_INT reflected_node=node;reflected_node[axis]=reflection_times_two-node[axis];
            u_ghost(node)=u_ghost(reflected_node);}}
}
//#####################################################################
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<float,1> >,float>;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<float,2> >,float>;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<float,3> >,float>;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<float,1> >,MATRIX<float,1> >;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<float,2> >,SYMMETRIC_MATRIX<float,2> >;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<float,3> >,SYMMETRIC_MATRIX<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<double,1> >,double>;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<double,2> >,double>;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<double,3> >,double>;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<double,1> >,MATRIX<double,1> >;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<double,2> >,SYMMETRIC_MATRIX<double,2> >;
template class BOUNDARY_REFLECTION_UNIFORM<GRID<VECTOR<double,3> >,SYMMETRIC_MATRIX<double,3> >;
#endif
