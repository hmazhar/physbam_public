#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <cstring>
using namespace PhysBAM;
//#####################################################################
template<class T> const int MAP_OCTREE_MESH<T>::
nodes_by_axis[3][2][4]={{{1,3,5,7},{0,2,4,6}},{{2,3,6,7},{0,1,4,5}},{{4,5,6,7},{0,1,2,3}}};
//#####################################################################
template<class T> const int MAP_OCTREE_MESH<T>::
face_by_axis[3][2]={{1,0},{3,2},{5,4}};
//#####################################################################
template<class T> const int MAP_OCTREE_MESH<T>::
edge_nodes_by_axis[3][4][2]={{{6,7},{2,3},{0,1},{4,5}},{{5,7},{4,6},{0,2},{1,3}},{{3,7},{1,5},{0,4},{2,6}}};
//#####################################################################
template<class T> const int MAP_OCTREE_MESH<T>::
opposite_face[6]={1,0,3,2,5,4};
//#####################################################################
// Function Map_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,
        n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    Map_Nodes(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
}
template<class T> void MAP_OCTREE_MESH<T>::
Map_Nodes(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    // process the nodes
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Node_Process_Node(pointer,cells(i,j,ij),cells(i+1,j,ij),cells(i,j+1,ij),cells(i+1,j+1,ij),cells(i,j,ij+1),cells(i+1,j,ij+1),cells(i,j+1,ij+1),cells(i+1,j+1,ij+1),map_function);
    // process the edges
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Node_Process_Edge(pointer,cells(i,j,ij),cells(i,j,ij+1),cells(i,j+1,ij+1),cells(i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Node_Process_Edge(pointer,cells(i,j,ij),cells(i+1,j,ij),cells(i+1,j,ij+1),cells(i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Node_Process_Edge(pointer,cells(i,j,ij),cells(i,j+1,ij),cells(i+1,j+1,ij),cells(i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    // process the faces
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Node_Process_Face(pointer,cells(i,j,ij),cells(i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Node_Process_Face(pointer,cells(i,j,ij),cells(i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Node_Process_Face(pointer,cells(i,j,ij),cells(i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    // process the cells
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Node_Process_Cell(pointer,cells(i,j,ij),map_function);
}
//#####################################################################
// Function Map_Edges
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Edges(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,
        n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    Map_Edges(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
}
template<class T> void MAP_OCTREE_MESH<T>::
Map_Edges(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis))
{
    // process the edges
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Edge_Process_Edge(pointer,cells(i,j,ij),cells(i,j,ij+1),cells(i,j+1,ij+1),cells(i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Edge_Process_Edge(pointer,cells(i,j,ij),cells(i+1,j,ij),cells(i+1,j,ij+1),cells(i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Edge_Process_Edge(pointer,cells(i,j,ij),cells(i,j+1,ij),cells(i+1,j+1,ij),cells(i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    // process the faces
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Edge_Process_Face(pointer,cells(i,j,ij),cells(i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Edge_Process_Face(pointer,cells(i,j,ij),cells(i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Edge_Process_Face(pointer,cells(i,j,ij),cells(i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    // process the cells
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Edge_Process_Cell(pointer,cells(i,j,ij),map_function);
}
//#####################################################################
// Function Map_Faces_For_Reduction
//#####################################################################
template<class T> template<class T3> void MAP_OCTREE_MESH<T>::
Map_Faces_For_Reduction(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,ARRAY<T3>& face_values,
          T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2,const T3 face_value_3,const T3 face_value_4))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,
        n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    // the faces between cells in the coarse uniform mesh
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Face_Process_Face_For_Reduction<T3>(pointer,cells(i,j,ij),cells(i+1,j,ij),OCTREE_CELL<T>::X_AXIS,face_values,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Face_Process_Face_For_Reduction<T3>(pointer,cells(i,j,ij),cells(i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,face_values,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Face_Process_Face_For_Reduction<T3>(pointer,cells(i,j,ij),cells(i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,face_values,map_function);
    // process the interior faces
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Face_Process_Cell_For_Reduction<T3>(pointer,cells(i,j,ij),face_values,map_function);
}
//#####################################################################
// Function Map_Boundary_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    // process the border faces
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++){
        Map_Node_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++){
        Map_Node_Process_Face(pointer,cells(i,0,ij),cells(i,1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++){
        Map_Node_Process_Face(pointer,cells(i,j,0),cells(i,j,1),OCTREE_CELL<T>::Z_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i,j,grid.numbers_of_cells.z+1),OCTREE_CELL<T>::Z_AXIS,map_function);}
    
    // process the border edges (note that the corner and edge nodes are only processed once (the loop start and end make sure of this)
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=0;j<=grid.numbers_of_cells.y;j++){ // front and back
        Map_Node_Process_Edge(pointer,cells(i,j,0),cells(i,j,1),cells(i,j+1,1),cells(i,j+1,0),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Node_Process_Edge(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i,j,grid.numbers_of_cells.z+1),cells(i,j+1,grid.numbers_of_cells.z+1),cells(i,j+1,grid.numbers_of_cells.z),OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z-1;ij++){ // bottom and top
        Map_Node_Process_Edge(pointer,cells(i,0,ij),cells(i,0,ij+1),cells(i,1,ij+1),cells(i,1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Node_Process_Edge(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y,ij+1),cells(i,grid.numbers_of_cells.y+1,ij+1),cells(i,grid.numbers_of_cells.y+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=0;ij<=grid.numbers_of_cells.z;ij++){ // left and right
        Map_Node_Process_Edge(pointer,cells(0,j,ij),cells(1,j,ij),cells(1,j,ij+1),cells(0,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
        Map_Node_Process_Edge(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),cells(grid.numbers_of_cells.x+1,j,ij+1),cells(grid.numbers_of_cells.x,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x-1;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++){ // front and back
        Map_Node_Process_Edge(pointer,cells(i,j,0),cells(i+1,j,0),cells(i+1,j,1),cells(i,j,1),OCTREE_CELL<T>::Y_AXIS,map_function);
        Map_Node_Process_Edge(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i+1,j,grid.numbers_of_cells.z),cells(i+1,j,grid.numbers_of_cells.z+1),cells(i,j,grid.numbers_of_cells.z+1),OCTREE_CELL<T>::Y_AXIS,map_function);}
    for(int j=0;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++){ // left and right
        Map_Node_Process_Edge(pointer,cells(0,j,ij),cells(0,j+1,ij),cells(1,j+1,ij),cells(1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
        Map_Node_Process_Edge(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x,j+1,ij),cells(grid.numbers_of_cells.x+1,j+1,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x-1;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++){ // bottom and top
        Map_Node_Process_Edge(pointer,cells(i,0,ij),cells(i,1,ij),cells(i+1,1,ij),cells(i+1,0,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
        Map_Node_Process_Edge(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),cells(i+1,grid.numbers_of_cells.y+1,ij),cells(i+1,grid.numbers_of_cells.y,ij),OCTREE_CELL<T>::Z_AXIS,map_function);}

    // process the border nodes (note that the corner and edge nodes are only processed once (the loop start and end make sure of this)
    for(int j=0;j<=grid.numbers_of_cells.y;j++)for(int ij=0;ij<=grid.numbers_of_cells.z;ij++){
        Map_Node_Process_Node(pointer,cells(0,j,ij),cells(1,j,ij),cells(0,j+1,ij),cells(1,j+1,ij),cells(0,j,ij+1),cells(1,j,ij+1),cells(0,j+1,ij+1),cells(1,j+1,ij+1),map_function);
        Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),cells(grid.numbers_of_cells.x,j+1,ij),cells(grid.numbers_of_cells.x+1,j+1,ij),cells(grid.numbers_of_cells.x,j,ij+1),cells(grid.numbers_of_cells.x+1,j,ij+1),cells(grid.numbers_of_cells.x,j+1,ij+1),cells(grid.numbers_of_cells.x+1,j+1,ij+1),map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x-1;i++)for(int ij=0;ij<=grid.numbers_of_cells.z;ij++){
        Map_Node_Process_Node(pointer,cells(i,0,ij),cells(i+1,0,ij),cells(i,1,ij),cells(i+1,1,ij),cells(i,0,ij+1),cells(i+1,0,ij+1),cells(i,1,ij+1),cells(i+1,1,ij+1),map_function);
        Map_Node_Process_Node(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i+1,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),cells(i+1,grid.numbers_of_cells.y+1,ij),cells(i,grid.numbers_of_cells.y,ij+1),cells(i+1,grid.numbers_of_cells.y,ij+1),cells(i,grid.numbers_of_cells.y+1,ij+1),cells(i+1,grid.numbers_of_cells.y+1,ij+1),map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x-1;i++)for(int j=1;j<=grid.numbers_of_cells.y-1;j++){
        Map_Node_Process_Node(pointer,cells(i,j,0),cells(i+1,j,0),cells(i,j+1,0),cells(i+1,j+1,0),cells(i,j,1),cells(i+1,j,1),cells(i,j+1,1),cells(i+1,j+1,1),map_function);
        Map_Node_Process_Node(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i+1,j,grid.numbers_of_cells.z),cells(i,j+1,grid.numbers_of_cells.z),cells(i+1,j+1,grid.numbers_of_cells.z),cells(i,j,grid.numbers_of_cells.z+1),cells(i+1,j,grid.numbers_of_cells.z+1),cells(i,j+1,grid.numbers_of_cells.z+1),cells(i+1,j+1,grid.numbers_of_cells.z+1),map_function);}
}
//#####################################################################
// Function Map_Left_Boundary_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Left_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                               void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int j=0;j<=grid.numbers_of_cells.y;j++)for(int ij=0;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Node(pointer,cells(0,j,ij),cells(1,j,ij),cells(0,j+1,ij),cells(1,j+1,ij),cells(0,j,ij+1),cells(1,j,ij+1),cells(0,j+1,ij+1),cells(1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Right_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                                void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int j=0;j<=grid.numbers_of_cells.y;j++)for(int ij=0;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),cells(grid.numbers_of_cells.x,j+1,ij),cells(grid.numbers_of_cells.x+1,j+1,ij),cells(grid.numbers_of_cells.x,j,ij+1),cells(grid.numbers_of_cells.x+1,j,ij+1),cells(grid.numbers_of_cells.x,j+1,ij+1),cells(grid.numbers_of_cells.x+1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Bottom_Boundary_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Bottom_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Face(pointer,cells(i,0,ij),cells(i,1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=0;i<=grid.numbers_of_cells.x;i++)for(int ij=0;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Node(pointer,cells(i,0,ij),cells(i+1,0,ij),cells(i,1,ij),cells(i+1,1,ij),cells(i,0,ij+1),cells(i+1,0,ij+1),cells(i,1,ij+1),cells(i+1,1,ij+1),map_function);
}
//#####################################################################
// Function Map_Top_Boundary_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Top_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Face(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=0;i<=grid.numbers_of_cells.x;i++)for(int ij=0;ij<=grid.numbers_of_cells.z;ij++)Map_Node_Process_Node(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i+1,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),cells(i+1,grid.numbers_of_cells.y+1,ij),cells(i,grid.numbers_of_cells.y,ij+1),cells(i+1,grid.numbers_of_cells.y,ij+1),cells(i,grid.numbers_of_cells.y+1,ij+1),cells(i+1,grid.numbers_of_cells.y+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Front_Boundary_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Front_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Face(pointer,cells(i,j,0),cells(i,j,1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=0;i<=grid.numbers_of_cells.x;i++)for(int j=0;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Node(pointer,cells(i,j,0),cells(i+1,j,0),cells(i,j+1,0),cells(i+1,j+1,0),cells(i,j,1),cells(i+1,j,1),cells(i,j+1,1),cells(i+1,j+1,1),map_function);
}
//#####################################################################
// Function Map_Back_Boundary_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Back_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Face(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i,j,grid.numbers_of_cells.z+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=0;i<=grid.numbers_of_cells.x;i++)for(int j=0;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Node(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i+1,j,grid.numbers_of_cells.z),cells(i,j+1,grid.numbers_of_cells.z),cells(i+1,j+1,grid.numbers_of_cells.z),cells(i,j,grid.numbers_of_cells.z+1),cells(i+1,j,grid.numbers_of_cells.z+1),cells(i,j+1,grid.numbers_of_cells.z+1),cells(i+1,j+1,grid.numbers_of_cells.z+1),map_function);
}
//#####################################################################
// Function Map_Ghost_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                     void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    Map_Left_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Right_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Bottom_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Top_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Front_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Back_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Ghost_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Left_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,
        m_end=-1,n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    for(int i=m_start;i<=m_end+1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Cell(pointer,Get_Cell(cells,i,j,ij),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end+1;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end+1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start;i<=m_end+1;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Right_Ghost_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Right_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                         void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int m_start=grid.numbers_of_cells.x+2,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,
        m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Cell(pointer,Get_Cell(cells,i,j,ij),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Bottom_Ghost_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Bottom_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,
        m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=-1,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end+1;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Cell(pointer,Get_Cell(cells,i,j,ij),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end+1;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end+1;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end+1;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Top_Ghost_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Top_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int m_start=1-number_of_ghost_cells,n_start=grid.numbers_of_cells.y+2,mn_start=1-number_of_ghost_cells,
        m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Cell(pointer,Get_Cell(cells,i,j,ij),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Front_Ghost_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Front_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                         void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,
        m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=-1;
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end+1;ij++) 
        Map_Node_Process_Cell(pointer,Get_Cell(cells,i,j,ij),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end+1;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end+1;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end+1;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Back_Ghost_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Back_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=grid.numbers_of_cells.z+2,
        m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Cell(pointer,Get_Cell(cells,i,j,ij),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i,j+1,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j,ij+1),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Edge(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i+1,j,ij),OCTREE_CELL<T>::Z_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)for(int ij=mn_start-1;ij<=mn_end;ij++) 
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
}
//#####################################################################
// Function Map_Exterior_Ghost_Cell_Nodes
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Exterior_Ghost_Cell_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                                                                     void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    int i,j,ij;
    // process the border faces
    for(j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Node_Process_Face(pointer,0,cells(1-number_of_ghost_cells,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells,j,ij),0,OCTREE_CELL<T>::X_AXIS,map_function);}
    for(i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Node_Process_Face(pointer,0,cells(i,1-number_of_ghost_cells,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(i,grid.numbers_of_cells.y+number_of_ghost_cells,ij),0,OCTREE_CELL<T>::Y_AXIS,map_function);}
    for(i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++){
        Map_Node_Process_Face(pointer,0,cells(i,j,1-number_of_ghost_cells),OCTREE_CELL<T>::Z_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(i,j,grid.numbers_of_cells.z+number_of_ghost_cells),0,OCTREE_CELL<T>::Z_AXIS,map_function);}
    // process the border nodes
    for(j=1-number_of_ghost_cells-1;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(ij=1-number_of_ghost_cells-1;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        i=1-number_of_ghost_cells-1;Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
        i=grid.numbers_of_cells.x+number_of_ghost_cells;Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);}
    for(i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells-1;i++)for(ij=1-number_of_ghost_cells-1;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        j=1-number_of_ghost_cells-1;Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
        j=grid.numbers_of_cells.y+number_of_ghost_cells;Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);}
    for(i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells-1;i++)for(j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells-1;j++){
        ij=1-number_of_ghost_cells-1;Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);
        ij=grid.numbers_of_cells.z+number_of_ghost_cells;Map_Node_Process_Node(pointer,Get_Cell(cells,i,j,ij),Get_Cell(cells,i+1,j,ij),Get_Cell(cells,i,j+1,ij),Get_Cell(cells,i+1,j+1,ij),Get_Cell(cells,i,j,ij+1),Get_Cell(cells,i+1,j,ij+1),Get_Cell(cells,i,j+1,ij+1),Get_Cell(cells,i+1,j+1,ij+1),map_function);}
}
//#####################################################################
// Function Map_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    Map_Left_Boundary_Faces(grid,cells,pointer,map_function);Map_Right_Boundary_Faces(grid,cells,pointer,map_function);
    Map_Bottom_Boundary_Faces(grid,cells,pointer,map_function);Map_Top_Boundary_Faces(grid,cells,pointer,map_function);
    Map_Front_Boundary_Faces(grid,cells,pointer,map_function);Map_Back_Boundary_Faces(grid,cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Left_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                               void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Face_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Right_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                                void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Bottom_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Bottom_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Face_Process_Face(pointer,cells(i,0,ij),cells(i,1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Top_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Top_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++)Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Front_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Front_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Face_Process_Face(pointer,cells(i,j,0),cells(i,j,1),OCTREE_CELL<T>::Z_AXIS,map_function);
}
//#####################################################################
// Function Map_Back_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Back_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Face_Process_Face(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i,j,grid.numbers_of_cells.z+1),OCTREE_CELL<T>::Z_AXIS,map_function);
}
//#####################################################################
// Function Map_Domain_Side_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Domain_Side_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int side,const int number_of_ghost_cells,void* pointer,
                                                               void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    switch(side){
      case 0:Map_Xmin_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 1:Map_Xmax_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 2:Map_Ymin_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 3:Map_Ymax_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 4:Map_Zmin_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 5:Map_Zmax_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;}
}
//#####################################################################
// Function Map_Xmin_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Xmin_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                                                               void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Face_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Xmax_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Xmax_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                                                                void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Ymin_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Ymin_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Face_Process_Face(pointer,cells(i,0,ij),cells(i,1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Top_Boundary_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Ymax_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Zmin_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Zmin_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)
        Map_Face_Process_Face(pointer,cells(i,j,0),cells(i,j,1),OCTREE_CELL<T>::Z_AXIS,map_function);
}
//#####################################################################
// Function Map_Zmax_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Zmax_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)
        Map_Face_Process_Face(pointer,cells(i,j,grid.numbers_of_cells.z),cells(i,j,grid.numbers_of_cells.z+1),OCTREE_CELL<T>::Z_AXIS,map_function);
}
//#####################################################################
// Function Map_Ghost_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Ghost_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                                                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    int m_start,n_start,mn_start,m_end,n_end,mn_end;
    n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    m_start=1-number_of_ghost_cells,m_end=0; Map_Faces(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
    m_start=grid.numbers_of_cells.x+1,m_end=grid.numbers_of_cells.x+number_of_ghost_cells; Map_Faces(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);

    m_start=1,mn_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    n_start=1-number_of_ghost_cells,n_end=0; Map_Faces(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
    n_start=grid.numbers_of_cells.y+1,n_end=grid.numbers_of_cells.y+number_of_ghost_cells; Map_Faces(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
    
    m_start=1,n_start=1,m_end=grid.numbers_of_cells.x,n_end=grid.numbers_of_cells.y;
    mn_start=1-number_of_ghost_cells,mn_end=0; Map_Faces(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
    mn_start=grid.numbers_of_cells.z+1,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells; Map_Faces(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);

    // map the gaps between the regions above
    for(int j=1-number_of_ghost_cells;j<=0;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Face_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int j=grid.numbers_of_cells.y+1;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Face_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1-number_of_ghost_cells;ij<=0;ij++){
        Map_Face_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=grid.numbers_of_cells.z+1;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Face_Process_Face(pointer,cells(0,j,ij),cells(1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j,ij),cells(grid.numbers_of_cells.x+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1-number_of_ghost_cells;ij<=0;ij++){
        Map_Face_Process_Face(pointer,cells(i,0,ij),cells(i,1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=grid.numbers_of_cells.z+1;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Face_Process_Face(pointer,cells(i,0,ij),cells(i,1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y,ij),cells(i,grid.numbers_of_cells.y+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);}
}
//#####################################################################
// Function Map_Exterior_Ghost_Cell_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Exterior_Ghost_Cell_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                                                                     void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    // process the border faces
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Face_Process_Face(pointer,0,cells(1-number_of_ghost_cells,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells,j,ij),0,OCTREE_CELL<T>::X_AXIS,map_function);}
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++){
        Map_Face_Process_Face(pointer,0,cells(i,1-number_of_ghost_cells,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y+number_of_ghost_cells,ij),0,OCTREE_CELL<T>::Y_AXIS,map_function);}
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++){
        Map_Face_Process_Face(pointer,0,cells(i,j,1-number_of_ghost_cells),OCTREE_CELL<T>::Z_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(i,j,grid.numbers_of_cells.z+number_of_ghost_cells),0,OCTREE_CELL<T>::Z_AXIS,map_function);}
}
//#####################################################################
// Function Map_Ghost_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Ghost_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
                                                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    for(int i=1-number_of_ghost_cells;i<=0;i++)for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Cell_Process_Cell(pointer,cells(i,j,ij),map_function);
    for(int i=grid.numbers_of_cells.x+1;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Cell_Process_Cell(pointer,cells(i,j,ij),map_function);

    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1-number_of_ghost_cells;j<=0;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Cell_Process_Cell(pointer,cells(i,j,ij),map_function);
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=grid.numbers_of_cells.y+1;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)for(int ij=1-number_of_ghost_cells;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Cell_Process_Cell(pointer,cells(i,j,ij),map_function);

    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1-number_of_ghost_cells;ij<=0;ij++)
        Map_Cell_Process_Cell(pointer,cells(i,j,ij),map_function);
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=grid.numbers_of_cells.z+1;ij<=grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
        Map_Cell_Process_Cell(pointer,cells(i,j,ij),map_function);
}
//#####################################################################
// Function Map_Boundary_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    Map_Left_Boundary_Cells(grid,cells,pointer,map_function);Map_Right_Boundary_Cells(grid,cells,pointer,map_function);
    Map_Bottom_Boundary_Cells(grid,cells,pointer,map_function);Map_Top_Boundary_Cells(grid,cells,pointer,map_function);
    Map_Front_Boundary_Cells(grid,cells,pointer,map_function);Map_Back_Boundary_Cells(grid,cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Boundary_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Left_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                               void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++) Map_Cell_Process_Cell(pointer,cells(0,j,ij),map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Right_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                                                                void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++) Map_Cell_Process_Cell(pointer,cells(grid.numbers_of_cells.x+1,j,ij),map_function);
}
//#####################################################################
// Function Map_Bottom_Boundary_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Bottom_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++) Map_Cell_Process_Cell(pointer,cells(i,0,ij),map_function);
}
//#####################################################################
// Function Map_Top_Boundary_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Top_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int ij=1;ij<=grid.numbers_of_cells.z;ij++) Map_Cell_Process_Cell(pointer,cells(i,grid.numbers_of_cells.y+1,ij),map_function);
}
//#####################################################################
// Function Map_Front_Boundary_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Front_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++) Map_Cell_Process_Cell(pointer,cells(i,j,0),map_function);
}
//#####################################################################
// Function Map_Back_Boundary_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Back_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
                       void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++) Map_Cell_Process_Cell(pointer,cells(i,j,grid.numbers_of_cells.z+1),map_function);
}
//#####################################################################
// Function Map_Node_Process_Cell
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Node_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    if(!cell->Has_Children()) return;
    // recur on each 8 of the children of the cell (if there are any)
    for(int i=0;i<8;i++) Map_Node_Process_Cell(pointer,cell->Child(i),map_function);
    // recur on each of the internal faces of the cell
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    // recur on each of the edges (right hand rule for ordering)
    Map_Node_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),
        OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),
        OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),
        OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Node_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),
        OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Node_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),
        OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Node_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),
        OCTREE_CELL<T>::Z_AXIS,map_function);
    // recur on the nodes
    Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::LX_LY_LZ,cell),Child(OCTREE_CELL<T>::HX_LY_LZ,cell),Child(OCTREE_CELL<T>::LX_HY_LZ,cell),Child(OCTREE_CELL<T>::HX_HY_LZ,cell),
        Child(OCTREE_CELL<T>::LX_LY_HZ,cell),Child(OCTREE_CELL<T>::HX_LY_HZ,cell),Child(OCTREE_CELL<T>::LX_HY_HZ,cell),Child(OCTREE_CELL<T>::HX_HY_HZ,cell),map_function);
}
//#####################################################################
// Function Map_Node_Process_Face
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Node_Process_Face(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    // if neither of the trees have children,then we can stop recurring
    if(!Has_Children(c1) && !Has_Children(c2)) return;
    // at least one of the tree has children,so recur on the faces
    switch(orient){
        case OCTREE_CELL<T>::X_AXIS:
            // recur on the faces
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            // recur on the edges
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_HZ,c1),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),
                Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),map_function);
            break;
        case OCTREE_CELL<T>::Y_AXIS:
            // recur on the faces
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            // recur on the edges
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_LZ,c1),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),
                Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),map_function);
            break;
        case OCTREE_CELL<T>::Z_AXIS:
            // recur on the faces
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            // recur on the edges
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_HZ,c1),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),
                Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),map_function);
            break;}
}
//#####################################################################
// Function Map_Node_Process_Edge
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Node_Process_Edge(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const int orient,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    // check if we have a minimal edge
    if(!Has_Children(c1) && !Has_Children(c2) && !Has_Children(c3) && !Has_Children(c4)) return;
    // do not have a minimal edge,so recur on the edges
    switch(orient){
        case OCTREE_CELL<T>::X_AXIS:
            // recur on the edges
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::LX_LY_HZ,c4),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c3),Child(OCTREE_CELL<T>::HX_LY_HZ,c4),
                OCTREE_CELL<T>::X_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c4),Child(OCTREE_CELL<T>::HX_LY_HZ,c4),
                Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::HX_LY_LZ,c3),map_function);
            break;
        case OCTREE_CELL<T>::Y_AXIS:
            // recur on the edges
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::HX_LY_LZ,c4),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),Child(OCTREE_CELL<T>::LX_HY_LZ,c3),Child(OCTREE_CELL<T>::HX_HY_LZ,c4),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),
                Child(OCTREE_CELL<T>::HX_LY_LZ,c4),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::HX_HY_LZ,c4),Child(OCTREE_CELL<T>::LX_HY_LZ,c3),map_function);
            break;
        case OCTREE_CELL<T>::Z_AXIS:
            // recur on the edges
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::LX_HY_LZ,c4),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Node_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_HZ,c3),Child(OCTREE_CELL<T>::LX_HY_HZ,c4),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c4),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),
                Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c4),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_HZ,c3),map_function);
            break;}
}
//#####################################################################
// Function Map_Node_Process_Node
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Node_Process_Node(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8))
{
    // check if we have a minimal node
    if(!Has_Children(c1) && !Has_Children(c2) && !Has_Children(c3) && !Has_Children(c4) && !Has_Children(c5) && !Has_Children(c6) && !Has_Children(c7) && !Has_Children(c8)){
        map_function(pointer,c1,c2,c3,c4,c5,c6,c7,c8);return;}
    // do not have a minimal node,so recur on the nodes
    Map_Node_Process_Node(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_HZ,c3),Child(OCTREE_CELL<T>::LX_LY_HZ,c4),
        Child(OCTREE_CELL<T>::HX_HY_LZ,c5),Child(OCTREE_CELL<T>::LX_HY_LZ,c6),Child(OCTREE_CELL<T>::HX_LY_LZ,c7),Child(OCTREE_CELL<T>::LX_LY_LZ,c8),map_function);
}
//#####################################################################
// Function Map_Edge_Process_Cell
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Edge_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis))
{
    if(!cell->Has_Children()) return;
    // recur on each 8 of the children of the cell (if there are any)
    for(int i=0;i<8;i++) Map_Edge_Process_Cell(pointer,cell->Child(i),map_function);
    // recur on each of the internal faces of the cell
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Edge_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    // recur on each of the edges (right hand rule for ordering
    Map_Edge_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),
        OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Edge_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),
        OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Edge_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),
        OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Edge_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),
        OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Edge_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),
        OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Edge_Process_Edge(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),
        OCTREE_CELL<T>::Z_AXIS,map_function);
}
//#####################################################################
// Function Map_Edge_Process_Face
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Edge_Process_Face(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis))
{
    // if neither of the trees have children,then we can stop recurring
    if(!Has_Children(c1) && !Has_Children(c2)) return;
    // at least one of the tree has children,so recur on the faces
    switch (orient){
        case OCTREE_CELL<T>::X_AXIS:
            // recur on the faces
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            // recur on the edges
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_HZ,c1),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            break;
        case OCTREE_CELL<T>::Y_AXIS:
            // recur on the faces
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            // recur on the edges
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_LZ,c1),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            break;
        case OCTREE_CELL<T>::Z_AXIS:
            // recur on the faces
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Edge_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            // recur on the edges
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_HZ,c1),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            break;}
}
//#####################################################################
// Function Map_Edge_Process_Edge
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Edge_Process_Edge(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const int orient,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis))
{
    // check if we have a minimal edge
    if(!Has_Children(c1) && !Has_Children(c2) && !Has_Children(c3) && !Has_Children(c4)){map_function(pointer,c1,c2,c3,c4,orient);return;}
    // do not have a minimal edge,so recur on the edges
    switch(orient){
        case OCTREE_CELL<T>::X_AXIS:
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::LX_LY_HZ,c4),
                OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),Child(OCTREE_CELL<T>::HX_LY_LZ,c3),Child(OCTREE_CELL<T>::HX_LY_HZ,c4),
                OCTREE_CELL<T>::X_AXIS,map_function);
            break;
        case OCTREE_CELL<T>::Y_AXIS:
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::HX_LY_LZ,c4),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),Child(OCTREE_CELL<T>::LX_HY_LZ,c3),Child(OCTREE_CELL<T>::HX_HY_LZ,c4),
                OCTREE_CELL<T>::Y_AXIS,map_function);
            break;
        case OCTREE_CELL<T>::Z_AXIS:
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),Child(OCTREE_CELL<T>::LX_LY_LZ,c3),Child(OCTREE_CELL<T>::LX_HY_LZ,c4),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Edge_Process_Edge(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),Child(OCTREE_CELL<T>::LX_LY_HZ,c3),Child(OCTREE_CELL<T>::LX_HY_HZ,c4),
                OCTREE_CELL<T>::Z_AXIS,map_function);
            break;}
}
//#####################################################################
// Function Map_Face_Process_Cell
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Face_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    if (!cell->Has_Children()) return;
    // recur on each 8 of the children of the cell (if there are any)
    for (int i=0;i<8;i++) Map_Face_Process_Cell(pointer,cell->Child(i),map_function);
    // recur on each of the internal faces of the cell
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::X_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,map_function);
}
//#####################################################################
// Function Map_Face_Process_Face
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Face_Process_Face(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    // if neither of the trees have children,then we can stop recurring
    if(!Has_Children(c1) && !Has_Children(c2)){map_function(pointer,c1,c2,orient);return;}
    // at least one of the tree has children,so recur on the faces
    switch(orient){
        case OCTREE_CELL<T>::X_AXIS:
            // recur on the faces
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),OCTREE_CELL<T>::X_AXIS,map_function);
            break;
        case OCTREE_CELL<T>::Y_AXIS:
            // recur on the faces
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,map_function);
            break;
        case OCTREE_CELL<T>::Z_AXIS:
            // recur on the faces
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,map_function);
            break;}
}
//#####################################################################
// Function Map_Face_Process_Cell
//#####################################################################
template<class T> template<class T3> void MAP_OCTREE_MESH<T>::
Map_Face_Process_Cell_For_Reduction(void* pointer,const OCTREE_CELL<T>* cell,ARRAY<T3>& face_values,
    T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2,const T3 face_value_3,const T3 face_value_4))
{
    if (!cell->Has_Children()) return;
    // recur on each 8 of the children of the cell (if there are any)
    for (int i=0;i<8;i++) Map_Face_Process_Cell_For_Reduction<T3>(pointer,cell->Child(i),face_values,map_function);
    // recur on each of the internal faces of the cell
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_LZ),OCTREE_CELL<T>::X_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::X_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::X_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::X_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_LZ),OCTREE_CELL<T>::Y_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_HZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_HZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Y_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_LY_LZ),cell->Child(OCTREE_CELL<T>::LX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::HX_LY_LZ),cell->Child(OCTREE_CELL<T>::HX_LY_HZ),OCTREE_CELL<T>::Z_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::LX_HY_LZ),cell->Child(OCTREE_CELL<T>::LX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(OCTREE_CELL<T>::HX_HY_LZ),cell->Child(OCTREE_CELL<T>::HX_HY_HZ),OCTREE_CELL<T>::Z_AXIS,face_values,map_function);
}
//#####################################################################
// Function Map_Face_Process_Face_For_Reduction
//#####################################################################
template<class T> template<class T3> T3 MAP_OCTREE_MESH<T>::
Map_Face_Process_Face_For_Reduction(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,ARRAY<T3>& face_values,
     T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2,const T3 face_value_3,const T3 face_value_4))
{
    // if neither of the trees have children,then we return the base value
    if(!Has_Children(c1) && !Has_Children(c2)){
        const OCTREE_CELL<T>* cells[]={c1,c2};int deepest_cell,deepest_depth;Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
        return face_values(cells[deepest_cell]->Face(face_by_axis[orient][deepest_cell]));}
    // at least one of the tree has children,so recur on the faces
    T3 value=0;
    switch(orient){
      case OCTREE_CELL<T>::X_AXIS:
        // recur on the faces
        value=map_function(pointer,Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_LY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::X_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::X_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::X_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_HZ,c2),OCTREE_CELL<T>::X_AXIS,face_values,map_function));
        break;
      case OCTREE_CELL<T>::Y_AXIS:
        // recur on the faces
        value=map_function(pointer,Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::LX_HY_LZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_HY_LZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Y_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_HZ,c2),OCTREE_CELL<T>::Y_AXIS,face_values,map_function));
        break;
      case OCTREE_CELL<T>::Z_AXIS:
        // recur on the faces
        value=map_function(pointer,Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::LX_LY_HZ,c1),Child(OCTREE_CELL<T>::LX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_LY_HZ,c1),Child(OCTREE_CELL<T>::HX_LY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::LX_HY_HZ,c1),Child(OCTREE_CELL<T>::LX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(OCTREE_CELL<T>::HX_HY_HZ,c1),Child(OCTREE_CELL<T>::HX_HY_LZ,c2),OCTREE_CELL<T>::Z_AXIS,face_values,map_function));
        break;}
    const OCTREE_CELL<T>* cells[]={c1,c2};int deepest_cell,deepest_depth;Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
    face_values(cells[deepest_cell]->Face(face_by_axis[orient][deepest_cell]))=value;
    return value;
}
//#####################################################################
// Function Map_Octree_Process_Cell
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Cell_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    if(!cell->Has_Children()) {map_function(pointer,cell);return;}
    for(int i=0;i<8;i++) Map_Cell_Process_Cell(pointer,cell->Child(i),map_function);
}
//#####################################################################
// Function Map_Cells
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,
        n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    Map_Cells(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
}
template<class T> void MAP_OCTREE_MESH<T>::
Map_Cells(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1))
{
    // process the interior cells
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Cell_Process_Cell(pointer,cells(i,j,ij),map_function);
}
//#####################################################################
// Function Map_Faces
//#####################################################################
template<class T> void MAP_OCTREE_MESH<T>::
Map_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,mn_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,
        n_end=grid.numbers_of_cells.y+number_of_ghost_cells,mn_end=grid.numbers_of_cells.z+number_of_ghost_cells;
    Map_Faces(m_start,m_end,n_start,n_end,mn_start,mn_end,cells,pointer,map_function);
}
template<class T> void MAP_OCTREE_MESH<T>::
Map_Faces(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
          void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis))
{
    // the faces between cells in the coarse uniform mesh
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Face_Process_Face(pointer,cells(i,j,ij),cells(i+1,j,ij),OCTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Face_Process_Face(pointer,cells(i,j,ij),cells(i,j+1,ij),OCTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end-1;ij++) Map_Face_Process_Face(pointer,cells(i,j,ij),cells(i,j,ij+1),OCTREE_CELL<T>::Z_AXIS,map_function);
    // process the interior faces
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) Map_Face_Process_Cell(pointer,cells(i,j,ij),map_function);
}
//#####################################################################
template class MAP_OCTREE_MESH<float>;
template void MAP_OCTREE_MESH<float>::Map_Faces_For_Reduction<float>(const GRID<VECTOR<float,3> >&,const ARRAY<OCTREE_CELL<float>*,VECTOR<int,3> >&,const int,void*,ARRAY<float>&, float (*map_function)(void*,const float,const float,const float,const float));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MAP_OCTREE_MESH<double>;
template void MAP_OCTREE_MESH<double>::Map_Faces_For_Reduction<double>(const GRID<VECTOR<double,3> >&,const ARRAY<OCTREE_CELL<double>*,VECTOR<int,3> >&,const int,void*,ARRAY<double>&, double (*map_function)(void*,const double,const double,const double,const double));
#endif
#endif
