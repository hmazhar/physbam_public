//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_UNIFORM  
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_UNIFORM<T_GRID,T2>::
BOUNDARY_UNIFORM()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_UNIFORM<T_GRID,T2>::
~BOUNDARY_UNIFORM()
{}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    T_ARRAYS_T2::Put(u,u_ghost);
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++)Fill_Single_Ghost_Region(grid,u_ghost,side,regions(side));
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid() && !clamp_below && !clamp_above);BOUNDARY_UNIFORM<T_GRID,T2> temp_boundary;
    if(use_fixed_boundary) temp_boundary.Set_Fixed_Boundary(true,fixed_boundary_value);
    for(int axis=1;axis<=T_GRID::dimension;axis++)temp_boundary.Fill_Ghost_Cells(grid.Get_Face_Grid(axis),u.Component(axis),u_ghost.Component(axis),0,time,number_of_ghost_cells);
}
//#####################################################################
// Function Find_Ghost_Regions
//#####################################################################
static inline void Find_Ghost_Regions_Helper(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& domain,const RANGE<VECTOR<int,1> >& ghost)
{
    regions(1)=RANGE<VECTOR<int,1> >(ghost.min_corner.x,domain.min_corner.x-1); // left
    regions(2)=RANGE<VECTOR<int,1> >(domain.max_corner.x+1,ghost.max_corner.x); // right
}
static inline void Find_Ghost_Regions_Helper(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& domain,const RANGE<VECTOR<int,2> >& ghost)
{
    regions(1)=RANGE<VECTOR<int,2> >(ghost.min_corner.x,domain.min_corner.x-1,domain.min_corner.y,domain.max_corner.y); // left
    regions(2)=RANGE<VECTOR<int,2> >(domain.max_corner.x+1,ghost.max_corner.x,domain.min_corner.y,domain.max_corner.y); // right
    regions(3)=RANGE<VECTOR<int,2> >(ghost.min_corner.x,ghost.max_corner.x,ghost.min_corner.y,domain.min_corner.y-1); // bottom
    regions(4)=RANGE<VECTOR<int,2> >(ghost.min_corner.x,ghost.max_corner.x,domain.max_corner.y+1,ghost.max_corner.y); // top
}
static inline void Find_Ghost_Regions_Helper(ARRAY<RANGE<VECTOR<int,3> > >& regions,const RANGE<VECTOR<int,3> >& domain,const RANGE<VECTOR<int,3> >& ghost)
{
    regions(1)=RANGE<VECTOR<int,3> >(ghost.min_corner.x,domain.min_corner.x-1,domain.min_corner.y,domain.max_corner.y,domain.min_corner.z,domain.max_corner.z); // left
    regions(2)=RANGE<VECTOR<int,3> >(domain.max_corner.x+1,ghost.max_corner.x,domain.min_corner.y,domain.max_corner.y,domain.min_corner.z,domain.max_corner.z); // right
    regions(3)=RANGE<VECTOR<int,3> >(ghost.min_corner.x,ghost.max_corner.x,ghost.min_corner.y,domain.min_corner.y-1,domain.min_corner.z,domain.max_corner.z); // bottom
    regions(4)=RANGE<VECTOR<int,3> >(ghost.min_corner.x,ghost.max_corner.x,domain.max_corner.y+1,ghost.max_corner.y,domain.min_corner.z,domain.max_corner.z); // top
    regions(5)=RANGE<VECTOR<int,3> >(ghost.min_corner.x,ghost.max_corner.x,ghost.min_corner.y,ghost.max_corner.y,ghost.min_corner.z,domain.min_corner.z-1); // front
    regions(6)=RANGE<VECTOR<int,3> >(ghost.min_corner.x,ghost.max_corner.x,ghost.min_corner.y,ghost.max_corner.y,domain.max_corner.z+1,ghost.max_corner.z); // back
}
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Find_Ghost_Regions(const T_GRID& grid,ARRAY<RANGE<TV_INT> >& regions,const int ghost_cells) const
{
    RANGE<TV_INT> domain=grid.Domain_Indices(),ghost=domain.Thickened(ghost_cells);
    regions.Resize(T_GRID::number_of_faces_per_cell);
    Find_Ghost_Regions_Helper(regions,domain,ghost);
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Fill_Single_Ghost_Region_Threaded(RANGE<TV_INT>& region,const T_GRID& grid,T_ARRAYS_T2& u_ghost,const int side)
{
    Fill_Single_Ghost_Region(grid,u_ghost,side,region);
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const int side,const RANGE<TV_INT>& region) const
{
    int axis=(side+1)/2,boundary=side&1?region.Maximum_Corner()[axis]+1:region.Minimum_Corner()[axis]-1;
    NODE_ITERATOR iterator(grid,region);
    if(use_fixed_boundary) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
        u_ghost(node)=fixed_boundary_value;}
    else if(clamp_below&&clamp_above) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=clamp(u_ghost(boundary_node),lower_threshold,upper_threshold);}
    else if(clamp_below) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=clamp_min(u_ghost(boundary_node),lower_threshold);}
    else if(clamp_above) for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=clamp_max(u_ghost(boundary_node),upper_threshold);}
    else for(;iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary;
        u_ghost(node)=u_ghost(boundary_node);}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time)
{
}
//#####################################################################
// Function Apply_Boundary_Condition_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_T2& u,const T time)
{
}
//#####################################################################
// Function Apply_Boundary_Condition_Single_Side
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Apply_Boundary_Condition_Single_Side(const T_GRID& grid,T_ARRAYS_T2& u,const int side,const T time) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_UNIFORM<T_GRID,T2>::
Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,1> >,float>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,1> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,2> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,3> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,4> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,2> >,float>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,2> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,3> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,4> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,5> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,3> >,float>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,3> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,5> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,6> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,3> >,bool>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,1> >,MATRIX<float,1> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,2> >,SYMMETRIC_MATRIX<float,2> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<float,3> >,SYMMETRIC_MATRIX<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,1> >,double>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,1> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,2> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,3> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,4> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,2> >,double>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,2> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,3> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,4> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,5> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,3> >,double>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,3> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,5> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,6> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,3> >,bool>;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,1> >,MATRIX<double,1> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,2> >,SYMMETRIC_MATRIX<double,2> >;
template class BOUNDARY_UNIFORM<GRID<VECTOR<double,3> >,SYMMETRIC_MATRIX<double,3> >;
#endif
