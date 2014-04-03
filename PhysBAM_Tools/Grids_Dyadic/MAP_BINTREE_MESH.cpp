//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <cstring>
using namespace PhysBAM;
//#####################################################################
template<class T> const int MAP_BINTREE_MESH<T>::
face_by_axis[2]={1,0};
//#####################################################################
template<class T> const int MAP_BINTREE_MESH<T>::
nodes_by_axis[][2][1]={{{0},{1}}};
//#####################################################################
template<class T> const int MAP_BINTREE_MESH<T>::
opposite_face[2]={1,0};
//#####################################################################
// Function Map_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    int m_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells;
    Map_Nodes(m_start,m_end,cells,pointer,map_function);
}
template<class T> void MAP_BINTREE_MESH<T>::
Map_Nodes(const int m_start,const int m_end,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer, void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    for(int i=m_start;i<=m_end-1;i++) Map_Node_Process_Node(pointer,cells(i),cells(i+1),map_function); // process the nodes
    for(int i=m_start;i<=m_end;i++)   Map_Node_Process_Cell(pointer,cells(i),map_function);            // process the cells
}
//#####################################################################
// Function Map_Faces_For_Reduction
//#####################################################################
template<class T> template<class T3> void MAP_BINTREE_MESH<T>::
Map_Faces_For_Reduction(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,ARRAY<T3>& face_values,T3 (*map_function)(void* pointer,const T3 face_value))
{
    int m_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells;
    for(int i=m_start;i<=m_end-1;i++) Map_Face_Process_Face_For_Reduction<T3>(pointer,cells(i),cells(i+1),face_values,map_function); // the faces between cells in the coarse uniform mesh
    for(int i=m_start;i<=m_end;i++)   Map_Face_Process_Cell_For_Reduction<T3>(pointer,cells(i),face_values,map_function);            // process the interior faces
}
//#####################################################################
// Function Map_Boundary_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                   void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    Map_Node_Process_Node(pointer,cells(0),                      cells(1),                        map_function);
    Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x),cells(grid.numbers_of_cells.x+1),map_function);
}
//#####################################################################
// Function Map_Left_Boundary_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Left_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    Map_Node_Process_Node(pointer,cells(0),                      cells(1),                        map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Right_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                         void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x),cells(grid.numbers_of_cells.x+1),map_function);
}
//#####################################################################
// Function Map_Ghost_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                     void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    Map_Left_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Right_Ghost_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
    Map_Exterior_Ghost_Cell_Nodes(grid,cells,number_of_ghost_cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Ghost_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Left_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                     void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    int m_start=1-number_of_ghost_cells,m_end=-1;
    for(int i=m_start;i<=m_end+1;i++) Map_Node_Process_Cell(pointer,cells(i),map_function);
    for(int i=m_start-1;i<=m_end;i++) Map_Node_Process_Node(pointer,Get_Cell(cells,i),Get_Cell(cells,i+1),map_function);
}
//#####################################################################
// Function Map_Right_Ghost_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Right_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    int m_start=grid.numbers_of_cells.x+2,m_end=grid.numbers_of_cells.x+number_of_ghost_cells;
    for(int i=m_start-1;i<=m_end;i++) Map_Node_Process_Cell(pointer,cells(i),map_function);
    for(int i=m_start-1;i<=m_end;i++) Map_Node_Process_Node(pointer,Get_Cell(cells,i),Get_Cell(cells,i+1),map_function);
}
//#####################################################################
// Function Map_Exterior_Ghost_Cell_Nodes
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Exterior_Ghost_Cell_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    Map_Node_Process_Node(pointer,0,                                                   cells(1-number_of_ghost_cells),map_function); // process the border nodes
    Map_Node_Process_Node(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells),0,                             map_function);
}
//#####################################################################
// Function Map_Boundary_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Boundary_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                   void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    Map_Left_Boundary_Faces(grid,cells,pointer,map_function);Map_Right_Boundary_Faces(grid,cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Boundary_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Left_Boundary_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    Map_Face_Process_Face(pointer,cells(0),cells(1),map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Right_Boundary_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x),cells(grid.numbers_of_cells.x+1),map_function);
}
//#####################################################################
// Function Map_Domain_Side_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Domain_Side_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int side,const int number_of_ghost_cells,void* pointer,
                                                               void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    if(side==0) Map_Xmin_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);
    else        Map_Xmax_Faces(grid,cells,number_of_ghost_cells,pointer,map_function);
}
//#####################################################################
// Function Map_Xmin_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Xmin_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    Map_Face_Process_Face(pointer,cells(0),cells(1),map_function);
}
//#####################################################################
// Function Map_Xmax_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Xmax_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x),cells(grid.numbers_of_cells.x+1),map_function);
}
//#####################################################################
// Function Map_Ghost_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Ghost_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                                                       void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    int m_start,m_end;
    m_start=1-number_of_ghost_cells,  m_end=0;                                             Map_Faces(m_start,m_end,cells,pointer,map_function);
    m_start=grid.numbers_of_cells.x+1,m_end=grid.numbers_of_cells.x+number_of_ghost_cells; Map_Faces(m_start,m_end,cells,pointer,map_function);
}
//#####################################################################
// Function Map_Exterior_Ghost_Cell_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Exterior_Ghost_Cell_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                              void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    // process the border faces
    Map_Face_Process_Face(pointer,0,                                                   cells(1-number_of_ghost_cells),map_function);
    Map_Face_Process_Face(pointer,cells(grid.numbers_of_cells.x+number_of_ghost_cells),0,                             map_function);
}
//#####################################################################
// Function Map_Ghost_Cells
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Ghost_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
                   void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1))
{
    for(int i=1-number_of_ghost_cells;i<=0;i++)
        Map_Cell_Process_Cell(pointer,cells(i),map_function);
    for(int i=grid.numbers_of_cells.x+1;i<=grid.numbers_of_cells.x+number_of_ghost_cells;i++)
        Map_Cell_Process_Cell(pointer,cells(i),map_function);
}
//#####################################################################
// Function Map_Boundary_Cells
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Boundary_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                   void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1))
{
    Map_Left_Boundary_Cells(grid,cells,pointer,map_function);Map_Right_Boundary_Cells(grid,cells,pointer,map_function);
}
//#####################################################################
// Function Map_Left_Boundary_Cells
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Left_Boundary_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1))
{
    Map_Cell_Process_Cell(pointer,cells(0),map_function);
}
//#####################################################################
// Function Map_Right_Boundary_Cells
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Right_Boundary_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
                        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1))
{
    Map_Cell_Process_Cell(pointer,cells(grid.numbers_of_cells.x+1),map_function);
}
//#####################################################################
// Function Map_Node_Process_Cell
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Node_Process_Cell(void* pointer,const BINTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    if(!cell->Has_Children()) return;
    // recur on each 2 of the children of the cell (if there are any)
    for(int i=0;i<2;i++) Map_Node_Process_Cell(pointer,cell->Child(i),map_function);
    // recur on the nodes
    Map_Node_Process_Node(pointer,Child(BINTREE_CELL<T>::LX,cell),Child(BINTREE_CELL<T>::HX,cell),map_function);
}
//#####################################################################
// Function Map_Node_Process_Node
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Node_Process_Node(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,
                      void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2))
{
    // check if we have a minimal node
    if(!Has_Children(c1) && !Has_Children(c2)){
        map_function(pointer,c1,c2);return;}
    // do not have a minimal node, so recur on the nodes
    Map_Node_Process_Node(pointer,Child(BINTREE_CELL<T>::HX,c1),Child(BINTREE_CELL<T>::LX,c2),map_function);
}
//#####################################################################
// Function Map_Face_Process_Cell
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Face_Process_Cell(void* pointer,const BINTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    if (!cell->Has_Children()) return;
    for (int i=0;i<2;i++) Map_Face_Process_Cell(pointer,cell->Child(i),map_function);
    Map_Face_Process_Face(pointer,cell->Child(BINTREE_CELL<T>::LX),cell->Child(BINTREE_CELL<T>::HX),map_function);
}
//#####################################################################
// Function Map_Face_Process_Cell_For_Reduction
//#####################################################################
template<class T> template<class T3> void MAP_BINTREE_MESH<T>::
Map_Face_Process_Cell_For_Reduction(void* pointer,const BINTREE_CELL<T>* cell,ARRAY<T3>& face_values,T3 (*map_function)(void* pointer,const T3 face_value))
{
    if (!cell->Has_Children()) return;
    // recur on both of the children of the cell (if there are any)
    for (int i=0;i<2;i++) Map_Face_Process_Cell_For_Reduction<T3>(pointer,cell->Child(i),face_values,map_function);
    // recur on each of the internal faces of the cell
    Map_Face_Process_Face_For_Reduction<T3>(pointer,cell->Child(BINTREE_CELL<T>::LX),cell->Child(BINTREE_CELL<T>::HX),face_values,map_function);
}
//#####################################################################
// Function Map_Face_Process_Face
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Face_Process_Face(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,
                      void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    if(!Has_Children(c1) && !Has_Children(c2)){map_function(pointer,c1,c2,0);return;}
    else Map_Face_Process_Face(pointer,Child(BINTREE_CELL<T>::HX,c1),Child(BINTREE_CELL<T>::LX,c2),map_function);
}
//#####################################################################
// Function Map_Face_Process_Face_For_Reduction
//#####################################################################
template<class T> template<class T3> T3 MAP_BINTREE_MESH<T>::
Map_Face_Process_Face_For_Reduction(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,ARRAY<T3>& face_values,
    T3 (*map_function)(void* pointer,const T3 face_value))
{
    // if neither of the trees have children,then we return the base value
    if(!Has_Children(c1) && !Has_Children(c2)){
        const BINTREE_CELL<T>* cells[]={c1,c2};int deepest_cell,deepest_depth;Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
        return face_values(cells[deepest_cell]->Face(face_by_axis[deepest_cell]));}
    // at least one of the tree has children,so recur on the faces
    T3 value=map_function(pointer,Map_Face_Process_Face_For_Reduction<T3>(pointer,Child(BINTREE_CELL<T>::HX,c1),Child(BINTREE_CELL<T>::LX,c2),face_values,map_function));

    const BINTREE_CELL<T>* cells[]={c1,c2};int deepest_cell,deepest_depth;Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
    face_values(cells[deepest_cell]->Face(face_by_axis[deepest_cell]))=value;
    return value;
}
//#####################################################################
// Function Map_Cell_Process_Cell
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Cell_Process_Cell(void* pointer,const BINTREE_CELL<T>* cell,
                      void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1))
{
    if(!cell->Has_Children()) {map_function(pointer,cell);return;}
    for(int i=0;i<2;i++) Map_Cell_Process_Cell(pointer,cell->Child(i),map_function);
}
//#####################################################################
// Function Map_Faces
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    int m_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells;
    Map_Faces(m_start,m_end,cells,pointer,map_function);
}
template<class T> void MAP_BINTREE_MESH<T>::
Map_Faces(const int m_start,const int m_end,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
          void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis))
{
    // the faces between cells in the coarse uniform mesh
    for(int i=m_start;i<=m_end-1;i++) Map_Face_Process_Face(pointer,cells(i),cells(i+1),map_function);
    // process the interior faces
    for(int i=m_start;i<=m_end;i++) Map_Face_Process_Cell(pointer,cells(i),map_function);
}
//#####################################################################
// Function Map_Cells
//#####################################################################
template<class T> void MAP_BINTREE_MESH<T>::
Map_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
          void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1))
{
    int m_start=1-number_of_ghost_cells,m_end=grid.numbers_of_cells.x+number_of_ghost_cells;
    Map_Cells(m_start,m_end,cells,pointer,map_function);
}
template<class T> void MAP_BINTREE_MESH<T>::
Map_Cells(const int m_start,const int m_end,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
          void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1))
{
    // process the interior cells
    for(int i=m_start;i<=m_end;i++) Map_Cell_Process_Cell(pointer,cells(i),map_function);
}
//#####################################################################
template class MAP_BINTREE_MESH<float>;
template void MAP_BINTREE_MESH<float>::Map_Faces_For_Reduction(const GRID<VECTOR<float,1> >&,const ARRAY<BINTREE_CELL<float>*,VECTOR<int,1> >&,const int,void*,ARRAY<float>&,float (*)(void*,const float));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MAP_BINTREE_MESH<double>;
template void MAP_BINTREE_MESH<double>::Map_Faces_For_Reduction(const GRID<VECTOR<double,1> >&,const ARRAY<BINTREE_CELL<double>*,VECTOR<int,1> >&,const int,void*,ARRAY<double>&,double (*)(void*,const double));
#endif

#endif
