//#####################################################################
// Copyright 2004-2008, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <cstring>
using namespace PhysBAM;
//#####################################################################
template<class T> const int MAP_QUADTREE_MESH<T>::
nodes_by_axis[][2][2]={{{1,3},{0,2}},{{2,3},{0,1}}};
//#####################################################################
template<class T> const int MAP_QUADTREE_MESH<T>::
face_by_axis[][2]={{1,0},{3,2}};
//#####################################################################
template<class T> const int MAP_QUADTREE_MESH<T>::
opposite_face[4]={1,0,3,2};
//#####################################################################
//#####################################################################
// Function Map_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    Map_Nodes(m_start,m_end,n_start,n_end,cells,pointer,map_function);
}
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Nodes(const int m_start,const int m_end,const int n_start,const int n_end,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
          void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    // process the nodes
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end-1;j++) Map_Node_Process_Node(pointer,cells(i,j),cells(i+1,j),cells(i,j+1),cells(i+1,j+1),map_function);
    // process the faces
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++) Map_Node_Process_Face(pointer,cells(i,j),cells(i+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++) Map_Node_Process_Face(pointer,cells(i,j),cells(i,j+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    // process the cells
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++) Map_Node_Process_Cell(pointer,cells(i,j),map_function);
}
//#####################################################################
// Function Map_Faces_For_Reduction
//#####################################################################
template<class T> template<class T3> void MAP_QUADTREE_MESH<T>::
Map_Faces_For_Reduction(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,ARRAY<T3>& face_values,T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    // the faces between cells in the coarse uniform mesh
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++) Map_Face_Process_Face_For_Reduction<T3>(pointer,cells(i,j),cells(i+1,j),QUADTREE_CELL<T>::X_AXIS,face_values,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++) Map_Face_Process_Face_For_Reduction<T3>(pointer,cells(i,j),cells(i,j+1),QUADTREE_CELL<T>::Y_AXIS,face_values,map_function);
    // process the interior faces
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++) Map_Face_Process_Cell_For_Reduction<T3>(pointer,cells(i,j),face_values,map_function);
}
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    // process the border faces
    for(int j=1;j<=grid.numbers_of_cells.y;j++){
        Map_Node_Process_Face(pointer,cells(0,j),cells(1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x;i++){
        Map_Node_Process_Face(pointer,cells(i,0),cells(i,1),QUADTREE_CELL<T>::Y_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(i,grid.numbers_of_cells.y),cells(i,grid.numbers_of_cells.y+1),QUADTREE_CELL<T>::Y_AXIS,map_function);}
    // process the border nodes
    for(int j=0;j<=grid.numbers_of_cells.y;j++){
        Map_Node_Process_Node(pointer,cells(0,j),cells(1,j),cells(0,j+1),cells(1,j+1),map_function);
        Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),cells(grid.numbers_of_cells.x,j+1),cells(grid.numbers_of_cells.x+1,j+1),map_function);}
    for(int i=1;i<=grid.numbers_of_cells.x-1;i++){
        Map_Node_Process_Node(pointer,cells(i,0),cells(i+1,0),cells(i,1),cells(i+1,1),map_function);
        Map_Node_Process_Node(pointer,cells(i,grid.numbers_of_cells.y),cells(i+1,grid.numbers_of_cells.y),cells(i,grid.numbers_of_cells.y+1),cells(i+1,grid.numbers_of_cells.y+1),map_function);}
}
//#####################################################################
// Function Map_Left_Boundary_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Left_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Face(pointer,cells(0,j),cells(1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int j=0;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Node(pointer,cells(0,j),cells(1,j),cells(0,j+1),cells(1,j+1),map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Right_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                         void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    // process the border faces
    for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Face(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int j=0;j<=grid.numbers_of_cells.y;j++)Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),cells(grid.numbers_of_cells.x,j+1),cells(grid.numbers_of_cells.x+1,j+1),map_function);
}
//#####################################################################
// Function Map_Bottom_Boundary_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Bottom_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                          void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)Map_Node_Process_Face(pointer,cells(i,0),cells(i,1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=0;i<=grid.numbers_of_cells.x;i++)Map_Node_Process_Node(pointer,cells(i,0),cells(i+1,0),cells(i,1),cells(i+1,1),map_function);
}
//#####################################################################
// Function Map_Top_Boundary_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Top_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                       void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)Map_Node_Process_Face(pointer,cells(i,grid.numbers_of_cells.y),cells(i,grid.numbers_of_cells.y+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=0;i<=grid.numbers_of_cells.x;i++)Map_Node_Process_Node(pointer,cells(i,grid.numbers_of_cells.y),cells(i+1,grid.numbers_of_cells.y),cells(i,grid.numbers_of_cells.y+1),cells(i+1,grid.numbers_of_cells.y+1),map_function);
}
//#####################################################################
// Function Map_Ghost_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                     void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    Map_Left_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Right_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Bottom_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Top_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Exterior_Ghost_Cell_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Ghost_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Left_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                     void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,m_end=-1,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    for(int i=m_start;i<=m_end+1;i++)for(int j=n_start;j<=n_end;j++)
        Map_Node_Process_Cell(pointer,cells(i,j),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end+1;i++)for(int j=n_start-1;j<=n_end;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i,j+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),Get_Cell(cells,i,j+1),Get_Cell(cells,i+1,j+1),map_function);
}
//#####################################################################
// Function Map_Right_Ghost_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Right_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    int m_start=grid.numbers_of_cells.x+2,n_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)
        Map_Node_Process_Cell(pointer,cells(i,j),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i,j+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),Get_Cell(cells,i,j+1),Get_Cell(cells,i+1,j+1),map_function);
}
//#####################################################################
// Function Map_Bottom_Ghost_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Bottom_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=-1;
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end+1;j++)
        Map_Node_Process_Cell(pointer,cells(i,j),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start;j<=n_end+1;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i,j+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),Get_Cell(cells,i,j+1),Get_Cell(cells,i+1,j+1),map_function);
}
//#####################################################################
// Function Map_Top_Ghost_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Top_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    int m_start=1-number_of_ghost_cells,n_start=grid.numbers_of_cells.y+2,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)
        Map_Node_Process_Cell(pointer,cells(i,j),map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++) 
        Map_Node_Process_Face(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i,j+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    for(int i=m_start-1;i<=m_end;i++)for(int j=n_start-1;j<=n_end;j++)
        Map_Node_Process_Node(pointer,Get_Cell(cells,i,j),Get_Cell(cells,i+1,j),Get_Cell(cells,i,j+1),Get_Cell(cells,i+1,j+1),map_function);
}
//#####################################################################
// Function Map_Exterior_Ghost_Cell_Nodes
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Exterior_Ghost_Cell_Nodes(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    // process the border faces
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++){
        Map_Node_Process_Face(pointer,0,cells(1-number_of_ghost_cells,j),QUADTREE_CELL<T>::X_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells,j),0,QUADTREE_CELL<T>::X_AXIS,map_function);}
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++){
        Map_Node_Process_Face(pointer,0,cells(i,1-number_of_ghost_cells),QUADTREE_CELL<T>::Y_AXIS,map_function);
        Map_Node_Process_Face(pointer,cells(i,grid.numbers_of_cells.y+number_of_ghost_cells),0,QUADTREE_CELL<T>::Y_AXIS,map_function);}
    // process the border nodes
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells-1;j++){
        Map_Node_Process_Node(pointer,0,cells(1-number_of_ghost_cells,j),0,cells(1-number_of_ghost_cells,j+1),map_function);
        Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells,j),0,cells(grid.numbers_of_cells.x+number_of_ghost_cells,j+1),0,map_function);}
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells-1;i++){
        Map_Node_Process_Node(pointer,0,0,cells(i,1-number_of_ghost_cells),cells(i+1,1-number_of_ghost_cells),map_function);
        Map_Node_Process_Node(pointer,cells(i,grid.numbers_of_cells.y+number_of_ghost_cells),cells(i+1,grid.numbers_of_cells.y+number_of_ghost_cells),0,0,map_function);}
    // corner nodes
    Map_Node_Process_Node(pointer,0,0,0,cells(1-number_of_ghost_cells,1-number_of_ghost_cells),map_function);
    Map_Node_Process_Node(pointer,0,0,cells(grid.numbers_of_cells.x+number_of_ghost_cells,1-number_of_ghost_cells),0,map_function);
    Map_Node_Process_Node(pointer,0,cells(1-number_of_ghost_cells,grid.numbers_of_cells.y+number_of_ghost_cells),0,0,map_function);
    Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells,grid.numbers_of_cells.y+number_of_ghost_cells),0,0,0,map_function);
}
//#####################################################################
// Function Map_Boundary_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Boundary_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    Map_Left_Boundary_Faces(grid,cells,pointer,map_function);Map_Right_Boundary_Faces(grid,cells,pointer,map_function);
    Map_Bottom_Boundary_Faces(grid,cells,pointer,map_function);Map_Top_Boundary_Faces(grid,cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Boundary_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Left_Boundary_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Face_Process_Face(pointer,cells(0,j),cells(1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Right_Boundary_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++)Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Bottom_Boundary_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Bottom_Boundary_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)Map_Face_Process_Face(pointer,cells(i,0),cells(i,1),QUADTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Top_Boundary_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Top_Boundary_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++)Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y),cells(i,grid.numbers_of_cells.y+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Domain_Side_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Domain_Side_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int side,const int number_of_ghost_cells,void* pointer,
                                                               void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    switch(side){
      case 0:Map_Xmin_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 1:Map_Xmax_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 2:Map_Ymin_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;
      case 3:Map_Ymax_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);break;}
}
//#####################################################################
// Function Map_Xmin_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Xmin_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)Map_Face_Process_Face(pointer,cells(0,j),cells(1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Xmax_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Xmax_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
}
//#####################################################################
// Function Map_Ymin_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Ymin_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)Map_Face_Process_Face(pointer,cells(i,0),cells(i,1),QUADTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Ymax_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Ymax_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y),cells(i,grid.numbers_of_cells.y+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Ghost_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Ghost_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                                                       void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    int m_start,n_start,m_end,n_end;
    n_start=1-number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    m_start=1-number_of_ghost_cells,m_end=0; Map_Faces(m_start,m_end,n_start,n_end,cells,pointer,map_function);
    m_start=grid.numbers_of_cells.x+1,m_end=grid.numbers_of_cells.x+number_of_ghost_cells; Map_Faces(m_start,m_end,n_start,n_end,cells,pointer,map_function);

    m_start=1,m_end=grid.numbers_of_cells.x;
    n_start=1-number_of_ghost_cells,n_end=0; Map_Faces(m_start,m_end,n_start,n_end,cells,pointer,map_function);
    n_start=grid.numbers_of_cells.y+1,n_end=grid.numbers_of_cells.y+number_of_ghost_cells; Map_Faces(m_start,m_end,n_start,n_end,cells,pointer,map_function);
    
    // map the gaps between the regions above
    for(int j=1-number_of_ghost_cells;j<=0;j++){
        Map_Face_Process_Face(pointer,cells(0,j),cells(1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);}
    for(int j=grid.numbers_of_cells.y+1;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++){
        Map_Face_Process_Face(pointer,cells(0,j),cells(1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x,j),cells(grid.numbers_of_cells.x+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);}
}
//#####################################################################
// Function Map_Exterior_Ghost_Cell_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Exterior_Ghost_Cell_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                              void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    // process the border faces
    for(int j=1-number_of_ghost_cells;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++){
        Map_Face_Process_Face(pointer,0,cells(1-number_of_ghost_cells,j),QUADTREE_CELL<T>::X_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells,j),0,QUADTREE_CELL<T>::X_AXIS,map_function);}
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++){
        Map_Face_Process_Face(pointer,0,cells(i,1-number_of_ghost_cells),QUADTREE_CELL<T>::Y_AXIS,map_function);
        Map_Face_Process_Face(pointer,cells(i,grid.numbers_of_cells.y+number_of_ghost_cells),0,QUADTREE_CELL<T>::Y_AXIS,map_function);}
}
//#####################################################################
// Function Map_Ghost_Cells
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Ghost_Cells(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    for(int i=1-number_of_ghost_cells;i<=0;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)
        Map_Cell_Process_Cell(pointer,cells(i,j),map_function);
    for(int i=grid.numbers_of_cells.x+1;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int j=1;j<=grid.numbers_of_cells.y;j++)
        Map_Cell_Process_Cell(pointer,cells(i,j),map_function);

    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int j=1-number_of_ghost_cells;j<=0;j++)
        Map_Cell_Process_Cell(pointer,cells(i,j),map_function);
    for(int i=1-number_of_ghost_cells;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)for(int j=grid.numbers_of_cells.y+1;j<=grid.numbers_of_cells.y+number_of_ghost_cells;j++)
        Map_Cell_Process_Cell(pointer,cells(i,j),map_function);
}
//#####################################################################
// Function Map_Boundary_Cells
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Boundary_Cells(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                   void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    Map_Left_Boundary_Cells(grid,cells,pointer,map_function);Map_Right_Boundary_Cells(grid,cells,pointer,map_function);
    Map_Bottom_Boundary_Cells(grid,cells,pointer,map_function);Map_Top_Boundary_Cells(grid,cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Boundary_Cells
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Left_Boundary_Cells(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++) Map_Cell_Process_Cell(pointer,cells(0,j),map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Cells
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Right_Boundary_Cells(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    for(int j=1;j<=grid.numbers_of_cells.y;j++) Map_Cell_Process_Cell(pointer,cells(grid.numbers_of_cells.x+1,j),map_function);
}
//#####################################################################
// Function Map_Bottom_Boundary_Cells
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Bottom_Boundary_Cells(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++) Map_Cell_Process_Cell(pointer,cells(i,0),map_function);
}
//#####################################################################
// Function Map_Top_Boundary_Cells
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Top_Boundary_Cells(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    for(int i=1;i<=grid.numbers_of_cells.x;i++) Map_Cell_Process_Cell(pointer,cells(i,grid.numbers_of_cells.y+1),map_function);
}
//#####################################################################
// Function Map_Node_Process_Cell
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Node_Process_Cell(void* pointer,const QUADTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    if(!cell->Has_Children()) return;
    // recur on each 4 of the children of the cell (if there are any)
    for(int i=0;i<4;i++) Map_Node_Process_Cell(pointer,cell->Child(i),map_function);
    // recur on each of the internal faces of the cell
    Map_Node_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::LX_LY),cell->Child(QUADTREE_CELL<T>::HX_LY),QUADTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::LX_HY),cell->Child(QUADTREE_CELL<T>::HX_HY),QUADTREE_CELL<T>::X_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::LX_LY),cell->Child(QUADTREE_CELL<T>::LX_HY),QUADTREE_CELL<T>::Y_AXIS,map_function);
    Map_Node_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::HX_LY),cell->Child(QUADTREE_CELL<T>::HX_HY),QUADTREE_CELL<T>::Y_AXIS,map_function);
    // recur on the nodes
    Map_Node_Process_Node(pointer,Child(QUADTREE_CELL<T>::LX_LY,cell),Child(QUADTREE_CELL<T>::HX_LY,cell),Child(QUADTREE_CELL<T>::LX_HY,cell),Child(QUADTREE_CELL<T>::HX_HY,cell),map_function);
}
//#####################################################################
// Function Map_Node_Process_Face
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Node_Process_Face(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const int orient,
                      void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    // if neither of the trees have children,then we can stop recurring
    if(!Has_Children(c1) && !Has_Children(c2)) return;
    // at least one of the tree has children, so recur on the faces
    switch(orient){
        case QUADTREE_CELL<T>::X_AXIS:
            // recur on the faces
            Map_Node_Process_Face(pointer,Child(QUADTREE_CELL<T>::HX_LY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),QUADTREE_CELL<T>::X_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::LX_HY,c2),QUADTREE_CELL<T>::X_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(QUADTREE_CELL<T>::HX_LY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),Child(QUADTREE_CELL<T>::HX_HY,c1),
                Child(QUADTREE_CELL<T>::LX_HY,c2),map_function);
            break;
        case QUADTREE_CELL<T>::Y_AXIS:
            // recur on the faces
            Map_Node_Process_Face(pointer,Child(QUADTREE_CELL<T>::LX_HY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),QUADTREE_CELL<T>::Y_AXIS,map_function);
            Map_Node_Process_Face(pointer,Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::HX_LY,c2),QUADTREE_CELL<T>::Y_AXIS,map_function);
            // recur on the nodes
            Map_Node_Process_Node(pointer,Child(QUADTREE_CELL<T>::LX_HY,c1),Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),
                Child(QUADTREE_CELL<T>::HX_LY,c2),map_function);
            break;}
}
//#####################################################################
// Function Map_Node_Process_Node
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Node_Process_Node(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4,
                      void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4))
{
    // check if we have a minimal node
    if(!Has_Children(c1) && !Has_Children(c2) && !Has_Children(c3) && !Has_Children(c4)){
        map_function(pointer,c1,c2,c3,c4);return;}
    // do not have a minimal node, so recur on the nodes
    Map_Node_Process_Node(pointer,Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::LX_HY,c2),Child(QUADTREE_CELL<T>::HX_LY,c3),
        Child(QUADTREE_CELL<T>::LX_LY,c4),map_function);
}
//#####################################################################
// Function Map_Face_Process_Cell
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Face_Process_Cell(void* pointer,const QUADTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    if (!cell->Has_Children()) return;
    // recur on each 4 of the children of the cell (if there are any)
    for (int i=0;i<4;i++) Map_Face_Process_Cell(pointer,cell->Child(i),map_function);
    // recur on each of the internal faces of the cell
    Map_Face_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::LX_LY),cell->Child(QUADTREE_CELL<T>::HX_LY),QUADTREE_CELL<T>::X_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::LX_HY),cell->Child(QUADTREE_CELL<T>::HX_HY),QUADTREE_CELL<T>::X_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::LX_LY),cell->Child(QUADTREE_CELL<T>::LX_HY),QUADTREE_CELL<T>::Y_AXIS,map_function);
    Map_Face_Process_Face(pointer,cell->Child(QUADTREE_CELL<T>::HX_LY),cell->Child(QUADTREE_CELL<T>::HX_HY),QUADTREE_CELL<T>::Y_AXIS,map_function);
}
//#####################################################################
// Function Map_Face_Process_Cell_For_Reduction
//#####################################################################
template<class T> template<class T3> void MAP_QUADTREE_MESH<T>::
Map_Face_Process_Cell_For_Reduction(void* pointer,const QUADTREE_CELL<T>* cell,ARRAY<T3>& face_values,T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2))
{
    if (!cell->Has_Children()) return;
    // recur on each 4 of the children of the cell (if there are any)
    for (int i=0;i<4;i++) Map_Face_Process_Cell_For_Reduction<T3>(pointer,cell->Child(i),face_values,map_function);
    // recur on each of the internal faces of the cell
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(QUADTREE_CELL<T>::LX_LY),cell->Child(QUADTREE_CELL<T>::HX_LY),QUADTREE_CELL<T>::X_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(QUADTREE_CELL<T>::LX_HY),cell->Child(QUADTREE_CELL<T>::HX_HY),QUADTREE_CELL<T>::X_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(QUADTREE_CELL<T>::LX_LY),cell->Child(QUADTREE_CELL<T>::LX_HY),QUADTREE_CELL<T>::Y_AXIS,face_values,map_function);
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(QUADTREE_CELL<T>::HX_LY),cell->Child(QUADTREE_CELL<T>::HX_HY),QUADTREE_CELL<T>::Y_AXIS,face_values,map_function);
}
//#####################################################################
// Function Map_Face_Process_Face
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Face_Process_Face(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const int orient,
                      void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    // if neither of the trees have children,then we can stop recurring
    if(!Has_Children(c1) && !Has_Children(c2)){map_function(pointer,c1,c2,orient);return;}
    // at least one of the tree has children,so recur on the faces
    switch(orient){
        case QUADTREE_CELL<T>::X_AXIS:
            // recur on the faces
            Map_Face_Process_Face(pointer,Child(QUADTREE_CELL<T>::HX_LY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),QUADTREE_CELL<T>::X_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::LX_HY,c2),QUADTREE_CELL<T>::X_AXIS,map_function);
            break;
        case QUADTREE_CELL<T>::Y_AXIS:
            // recur on the faces
            Map_Face_Process_Face(pointer,Child(QUADTREE_CELL<T>::LX_HY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),QUADTREE_CELL<T>::Y_AXIS,map_function);
            Map_Face_Process_Face(pointer,Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::HX_LY,c2),QUADTREE_CELL<T>::Y_AXIS,map_function);
            break;}
}
//#####################################################################
// Function Map_Face_Process_Face_For_Reduction
//#####################################################################
template<class T> template<class T3> T3 MAP_QUADTREE_MESH<T>::
Map_Face_Process_Face_For_Reduction(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const int orient,ARRAY<T3>& face_values,
    T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2))
{
    // if neither of the trees have children,then we return the base value
    if(!Has_Children(c1) && !Has_Children(c2)){
        const QUADTREE_CELL<T>* cells[]={c1,c2};int deepest_cell,deepest_depth;Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
        return face_values(cells[deepest_cell]->Face(face_by_axis[orient][deepest_cell]));}
    // at least one of the tree has children,so recur on the faces
    T3 value=0;
    switch(orient){
      case QUADTREE_CELL<T>::X_AXIS:
          // recur on the faces
        value=map_function(pointer,Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(QUADTREE_CELL<T>::HX_LY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),QUADTREE_CELL<T>::X_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::LX_HY,c2),QUADTREE_CELL<T>::X_AXIS,face_values,map_function));
        break;
      case QUADTREE_CELL<T>::Y_AXIS:
        // recur on the faces
        value=map_function(pointer,Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(QUADTREE_CELL<T>::LX_HY,c1),Child(QUADTREE_CELL<T>::LX_LY,c2),QUADTREE_CELL<T>::Y_AXIS,face_values,map_function),
                           Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(QUADTREE_CELL<T>::HX_HY,c1),Child(QUADTREE_CELL<T>::HX_LY,c2),QUADTREE_CELL<T>::Y_AXIS,face_values,map_function));
        break;}
    const QUADTREE_CELL<T>* cells[]={c1,c2};int deepest_cell,deepest_depth;Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
    face_values(cells[deepest_cell]->Face(face_by_axis[orient][deepest_cell]))=value;
    return value;
}
//#####################################################################
// Function Map_Cell_Process_Cell
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Cell_Process_Cell(void* pointer,const QUADTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    if(!cell->Has_Children()) {map_function(pointer,cell);return;}
    for(int i=0;i<4;i++) Map_Cell_Process_Cell(pointer,cell->Child(i),map_function);
}
//#####################################################################
// Function Map_Faces
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Faces(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    Map_Faces(m_start,m_end,n_start,n_end,cells,pointer,map_function);
}
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Faces(const int m_start,const int m_end,const int n_start,const int n_end,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
          void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,int axis))
{
    // the faces between cells in the coarse uniform mesh
    for(int i=m_start;i<=m_end-1;i++)for(int j=n_start;j<=n_end;j++) Map_Face_Process_Face(pointer,cells(i,j),cells(i+1,j),QUADTREE_CELL<T>::X_AXIS,map_function);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end-1;j++) Map_Face_Process_Face(pointer,cells(i,j),cells(i,j+1),QUADTREE_CELL<T>::Y_AXIS,map_function);
    // process the interior faces
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++) Map_Face_Process_Cell(pointer,cells(i,j),map_function);
}
//#####################################################################
// Function Map_Cells
//#####################################################################
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Cells(const GRID<TV>& grid,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    int m_start=1-number_of_ghost_cells,n_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells,n_end=grid.numbers_of_cells.y+number_of_ghost_cells;
    Map_Cells(m_start,m_end,n_start,n_end,cells,pointer,map_function);
}
template<class T> void MAP_QUADTREE_MESH<T>::
Map_Cells(const int m_start,const int m_end,const int n_start,const int n_end,const ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >& cells,void* pointer,
          void (*map_function)(void* pointer,const QUADTREE_CELL<T>* c1))
{
    // process the interior cells
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++) Map_Cell_Process_Cell(pointer,cells(i,j),map_function);
}
//#####################################################################
template class MAP_QUADTREE_MESH<float>;
template void MAP_QUADTREE_MESH<float>::Map_Faces_For_Reduction(const GRID<VECTOR<float,2> >&,const ARRAY<QUADTREE_CELL<float>*,VECTOR<int,2> >&,const int,void*,ARRAY<float>&,float (*)(void*,const float,const float));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MAP_QUADTREE_MESH<double>;
template void MAP_QUADTREE_MESH<double>::Map_Faces_For_Reduction(const GRID<VECTOR<double,2> >&,const ARRAY<QUADTREE_CELL<double>*,VECTOR<int,2> >&,const int,void*,ARRAY<double>&,double (*)(void*,const double,const double));
#endif

#endif
