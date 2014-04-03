//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_UNIFORM_PERIODIC
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM_PERIODIC.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM_PERIODIC<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    T_ARRAYS_T2::Put(u,u_ghost); // interior
    TV_INT periods=grid.Domain_Indices().Maximum_Corner()-TV_INT::All_Ones_Vector();
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){
        int side=2*axis+axis_side-2;
        TV_INT period=(axis_side==1?1:-1)*periods[axis]*TV_INT::Axis_Vector(axis);
        for(NODE_ITERATOR iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            u_ghost(node)=u_ghost(node+period);}}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM_PERIODIC<T_GRID,T2>::
Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time)
{
    assert(!grid.Is_MAC_Grid());
    for(int axis=1;axis<=T_GRID::dimension;axis++)
        for(NODE_ITERATOR iterator(grid,0,T_GRID::BOUNDARY_REGION,2*axis);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            TV_INT opposite_node=node;opposite_node[axis]=1;
            u(node)=u(opposite_node);}
}
