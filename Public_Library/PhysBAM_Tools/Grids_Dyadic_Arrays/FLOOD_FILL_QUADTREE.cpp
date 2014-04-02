//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_QUADTREE
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/FLOOD_FILL_QUADTREE.h>
using namespace PhysBAM;
//#####################################################################
// Function Flood_Fill
//#####################################################################
template<class T> static void Create_Cell_Index_To_Pointer_Array(void* data,const QUADTREE_CELL<T>* cell)
{
    ARRAY<const QUADTREE_CELL<T>*>* cell_index_to_pointer=(ARRAY<const QUADTREE_CELL<T>*>*)data;
    (*cell_index_to_pointer)(cell->Cell())=cell;
}
// colors should be initialized by the user with 0's where colors will be filled and -1 for uncolorable nodes (where no color will be filled)
// returns number of colors
template<class T> int FLOOD_FILL_QUADTREE<T>::
Flood_Fill(const QUADTREE_GRID<T>& grid,ARRAY<int>& colors,const ARRAY<bool>& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node)
{
    const QUADTREE_CELL<T>* seed_node;int fill_color=0; 
    cell_index_to_pointer.Resize(grid.number_of_cells);
    MAP_QUADTREE_MESH<T>::Map_Cells(grid.uniform_grid,grid.cells,grid.number_of_ghost_cells,&cell_index_to_pointer,Create_Cell_Index_To_Pointer_Array);
    flood_fill_stack.Preallocate(grid.number_of_cells);
    while(Find_Uncolored_Node(colors,seed_node)){
        bool touches_uncolorable_node;fill_color++;Flood_Fill_From_Seed_Node(grid,colors,fill_color,edge_is_blocked,touches_uncolorable_node,seed_node);
        if(color_touches_uncolorable_node)color_touches_uncolorable_node->Append(touches_uncolorable_node);}
    flood_fill_stack.Clean_Memory();
    return fill_color;
}
//#####################################################################
// Function Find_Uncolored_Node
//#####################################################################
template<class T> bool FLOOD_FILL_QUADTREE<T>::
Find_Uncolored_Node(const ARRAY<int>& colors,const QUADTREE_CELL<T>*& seed_cell)
{
    for(int i=1;i<=cell_index_to_pointer.m;i++)if(cell_index_to_pointer(i)&&colors(i)==0){seed_cell=cell_index_to_pointer(i);return true;}
    return false;
}
//#####################################################################
// Function Flood_Fill_From_Seed_Node
//#####################################################################
template<class T> void FLOOD_FILL_QUADTREE<T>::
Flood_Fill_From_Seed_Node(const QUADTREE_GRID<T>& grid,ARRAY<int>& colors,const int fill_color,const ARRAY<bool>& edge_is_blocked,bool& touches_uncolorable_node,const QUADTREE_CELL<T>* seed_node)
{
    assert(colors(seed_node->Cell())==0);touches_uncolorable_node=false;
    flood_fill_stack.Remove_All();
    flood_fill_stack.Push(seed_node);
    static const int opposite_face[]={1,0,3,2}; // gives the face opposite to the index
    while(!flood_fill_stack.Empty()){
        const QUADTREE_CELL<T>* node=flood_fill_stack.Pop();
        if(colors(node->Cell())==-1){touches_uncolorable_node=true;continue;}else if(colors(node->Cell())!=0)continue;colors(node->Cell())=fill_color;
        for(int i=0;i<4;i++){
            ARRAY<QUADTREE_CELL<T>*> face_neighbors;
            node->Get_All_Face_Neighbors(i,face_neighbors,&grid);
            for(int neighbor=1;neighbor<=face_neighbors.m;neighbor++){QUADTREE_CELL<T>* neighbor_node=face_neighbors(neighbor);
                bool blocked;if(node->Depth_Of_This_Cell()<=neighbor_node->Depth_Of_This_Cell())blocked=edge_is_blocked(neighbor_node->Face(opposite_face[i]));
                else blocked=edge_is_blocked(node->Face(i));
                if(!blocked&&colors(neighbor_node->Cell())<=0) flood_fill_stack.Push(neighbor_node);}}}
}
//#####################################################################
template class FLOOD_FILL_QUADTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FLOOD_FILL_QUADTREE<double>;
#endif
#endif
