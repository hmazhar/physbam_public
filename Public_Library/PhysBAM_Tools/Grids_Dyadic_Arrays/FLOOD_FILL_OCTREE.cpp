#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_OCTREE
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/FLOOD_FILL_OCTREE.h>
using namespace PhysBAM;
//#####################################################################
// Function Flood_Fill
//#####################################################################
// colors should be initialized by the user with 0's where colors will be filled and -1 for uncolorable nodes (where no color will be filled)
// returns number of colors
template<class T> int FLOOD_FILL_OCTREE<T>::
Flood_Fill(const OCTREE_GRID<T>& grid,ARRAY<int>& colors,const ARRAY<bool>& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node)
{
    const OCTREE_CELL<T>* seed_node;int fill_color=0;
    ARRAY<OCTREE_CELL<T>*>& cell_pointer_from_index=grid.Cell_Pointer_From_Index();
    flood_fill_stack.Preallocate(grid.number_of_cells);
    while(Find_Uncolored_Node(cell_pointer_from_index,colors,seed_node)){
        bool touches_uncolorable_node;fill_color++;Flood_Fill_From_Seed_Node(grid,colors,fill_color,edge_is_blocked,touches_uncolorable_node,seed_node);
        if(color_touches_uncolorable_node)color_touches_uncolorable_node->Append(touches_uncolorable_node);}
    flood_fill_stack.Clean_Memory();
    return fill_color;
}
//#####################################################################
// Function Find_Uncolored_Node
//#####################################################################
template<class T> bool FLOOD_FILL_OCTREE<T>::
Find_Uncolored_Node(const ARRAY<OCTREE_CELL<T>*>& cell_pointer_from_index,const ARRAY<int>& colors,const OCTREE_CELL<T>*& seed_cell)
{
    for(int i=1;i<=cell_pointer_from_index.m;i++)if(cell_pointer_from_index(i)&&colors(i)==0){seed_cell=cell_pointer_from_index(i);return true;}
    return false;
}
//#####################################################################
// Function Flood_Fill_From_Seed_Node
//#####################################################################
template<class T> void FLOOD_FILL_OCTREE<T>::
Flood_Fill_From_Seed_Node(const OCTREE_GRID<T>& grid,ARRAY<int>& colors,const int fill_color,const ARRAY<bool>& edge_is_blocked,bool& touches_uncolorable_node,const OCTREE_CELL<T>* seed_node)
{
    assert(colors(seed_node->Cell())==0);touches_uncolorable_node=false;
    flood_fill_stack.Remove_All();
    flood_fill_stack.Push(seed_node);
    static const int opposite_face[]={1,0,3,2,5,4}; // gives the face opposite to the index
    while(!flood_fill_stack.Empty()){
        const OCTREE_CELL<T>* node=flood_fill_stack.Pop();
        if(colors(node->Cell())==-1){touches_uncolorable_node=true;continue;}else if(colors(node->Cell())!=0)continue;colors(node->Cell())=fill_color;
        for(int i=0;i<6;i++){
            ARRAY<OCTREE_CELL<T>*> face_neighbors;
            node->Get_All_Face_Neighbors(i,face_neighbors,&grid);
            for(int neighbor=1;neighbor<=face_neighbors.m;neighbor++){OCTREE_CELL<T>* neighbor_node=face_neighbors(neighbor);
                bool blocked;if(node->Depth_Of_This_Cell()<=neighbor_node->Depth_Of_This_Cell())blocked=edge_is_blocked(neighbor_node->Face(opposite_face[i]));
                else blocked=edge_is_blocked(node->Face(i));
                if(!blocked&&colors(neighbor_node->Cell())<=0) flood_fill_stack.Push(neighbor_node);}}}
}
//#####################################################################
template class FLOOD_FILL_OCTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FLOOD_FILL_OCTREE<double>;
#endif
#endif
