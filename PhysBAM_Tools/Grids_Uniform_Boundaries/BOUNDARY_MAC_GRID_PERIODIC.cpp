//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_PERIODIC
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MAC_GRID_PERIODIC<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    T_ARRAYS_T2::Put(u,u_ghost); // interior
    TV_INT periods=grid.Domain_Indices().Maximum_Corner();
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){
        int side=2*axis+axis_side-2;
        TV_INT period=(axis_side==1?1:-1)*periods[axis]*TV_INT::Axis_Vector(axis);
        for(NODE_ITERATOR iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            u_ghost(node)=u_ghost(node+period);}}
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MAC_GRID_PERIODIC<T_GRID,T2>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid());
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        //Fill_Ghost_Cells(grid.Get_Face_Grid(axis,u.Component(axis),u_ghost.Component(axis),0,time,number_of_ghost_cells);
        const T_ARRAYS_T2& u_axis=u.Component(axis);
        T_ARRAYS_T2& u_ghost_axis=u_ghost.Component(axis);
        T_ARRAYS_T2::Put(u_axis,u_ghost_axis); // interior
        T_GRID face_grid=grid.Get_Face_Grid(axis);
        TV_INT periods=grid.Domain_Indices().Maximum_Corner();
        ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int face_axis=1;face_axis<=T_GRID::dimension;face_axis++)for(int axis_side=1;axis_side<=2;axis_side++){
            int side=2*face_axis+axis_side-2;
            TV_INT period=(axis_side==1?1:-1)*periods[face_axis]*TV_INT::Axis_Vector(face_axis);
            for(NODE_ITERATOR iterator(face_grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
                u_ghost_axis(node)=u_ghost_axis(node+period);}}}
}
//#####################################################################
template class BOUNDARY_MAC_GRID_PERIODIC<GRID<VECTOR<float,2> >,VECTOR<float,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<GRID<VECTOR<float,2> >,float>;
template class BOUNDARY_MAC_GRID_PERIODIC<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_MAC_GRID_PERIODIC<GRID<VECTOR<double,2> >,VECTOR<double,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<GRID<VECTOR<double,2> >,double>;
template class BOUNDARY_MAC_GRID_PERIODIC<GRID<VECTOR<double,3> >,double>;
#endif
