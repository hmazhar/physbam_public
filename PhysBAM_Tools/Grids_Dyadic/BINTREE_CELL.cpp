//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Cell_Compaction_Mapping
//#####################################################################
template<class T> void BINTREE_CELL<T>::
Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count)
{
    mapping_array(Cell())=++cell_count;Cell()=cell_count;
    if(Has_Children()) for(int i=0;i<2;i++) Child(i)->Create_Cell_Compaction_Mapping(mapping_array,cell_count);
}
//#####################################################################
// Function Storage_Requirement
//#####################################################################
template<class T> int BINTREE_CELL<T>::
Storage_Requirement() const
{
    int total=0;
    if(Has_Children()){total+=sizeof(BINTREE_CHILDREN<T>);for(int i=0;i<2;i++) total+=Child(i)->Storage_Requirement();}
    return total;   
}
//#####################################################################
// Function Create_Root
//#####################################################################
template<class T> BINTREE_CELL<T>* BINTREE_CELL<T>::
Create_Root(const bool use_nodes,const bool use_faces,int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,1>& root_center,const VECTOR<T,1>& root_DX)
{
    BINTREE_CHILDREN<T>* owner=new BINTREE_CHILDREN<T>(use_nodes,use_faces);
    owner->Initialize_Root_Cell(number_of_cells,number_of_nodes,number_of_faces,root_center,root_DX);
    return &(owner->children[0]);
}
//#####################################################################
// Function Create_Children
//#####################################################################
template<class T> void BINTREE_CELL<T>::
Create_Children(int& number_of_cells,ARRAY<BINTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_nodes,int& number_of_faces,
                ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_faces,const BINTREE_GRID<T>* grid)
{
    assert(children == 0);children=new BINTREE_CHILDREN<T>(owner->nodes!=0,owner->faces!=0);
    children->Initialize_Non_Root_Cell(this,number_of_cells,new_cells,number_of_nodes,new_nodes,number_of_faces,new_faces,Center(),grid);
}
//#####################################################################
// Function Delete_Children
//#####################################################################
// returns an array of the nodes that were deleted
template<class T> void BINTREE_CELL<T>::
Delete_Children()
{
    if(!Has_Children()) return;
    for(int i=0;i<2;i++) Child(i)->Delete_Children();
    delete children;children=0;
}
//#####################################################################
// Function Get_Neighbor
//#####################################################################
// returns the requested neighbor or 0 if no neighbor exists - offsets are between -1 and 1
// same_depth_only=true makes the function return a cell at the same depth or 0 - same_depth_only=false may return a neighbor cell higher up
template<class T> BINTREE_CELL<T>* BINTREE_CELL<T>::
Get_Neighbor(const int x_offset,const BINTREE_GRID<T>* grid,const ARRAY<VECTOR<BINTREE_CELL<T>*,2> >* neighbors,const bool same_depth_only) const
{
    assert(-1 == x_offset || x_offset == 1);
    if(neighbors) return (*neighbors)(Cell())(((x_offset+1)>>1)+1);
    if(owner->parent == 0){ // we are at the top of the tree... use the uniform grid
        if(grid==0) return 0;TV_INT i=grid->uniform_grid.Cell(Center(),grid->number_of_ghost_cells);
        if(grid->cells.Domain_Indices().Lazy_Outside(i+TV_INT(x_offset))) return 0;
        return grid->cells(i+TV_INT(x_offset));}
    int possible_child=child_index+x_offset;
    if(Can_Go_In_X_Direction(child_index,x_offset)){
        assert(possible_child >= 0 && possible_child < 2);return &owner->children[possible_child];}
    else{ // we need to keep going up the tree to find the node where we can go down the right branch
        BINTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(x_offset,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
        if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
        while(correct_branch->Has_Children()) correct_branch=correct_branch->Child(child_index-x_offset);
        return correct_branch;}
}
//#####################################################################
// Function Get_All_Face_Neighbors
//#####################################################################
template<class T> struct GET_ALL_FACE_NEIGHBORS_HELPER{int cell_to_add;ARRAY<BINTREE_CELL<T>*>* face_neighbors;};
template<class T> static void
Get_All_Face_Neighbors_Helper(void* data,const BINTREE_CELL<T>* cell1,const BINTREE_CELL<T>* cell2,const int axis)
{
    GET_ALL_FACE_NEIGHBORS_HELPER<T>* helper=(GET_ALL_FACE_NEIGHBORS_HELPER<T>*)data;
    const BINTREE_CELL<T>* cells[]={cell1,cell2};
    helper->face_neighbors->Append((BINTREE_CELL<T>*)cells[helper->cell_to_add]);
}
//#####################################################################
// Function Get_All_Face_Neighbors
//#####################################################################
template<class T> void BINTREE_CELL<T>::
Get_All_Face_Neighbors(int face_index,ARRAY<BINTREE_CELL<T>*>& face_neighbors,const BINTREE_GRID<T>* grid) const
{
    int face_index_to_offset[]={-1,1}; // convertes from face index to x offsets for use in neighbor search
    BINTREE_CELL<T>* neighbor=Get_Neighbor(face_index_to_offset[face_index],grid,0,false);
    if(!neighbor) return;
    if(neighbor->Has_Children()){
        GET_ALL_FACE_NEIGHBORS_HELPER<T> helper;helper.face_neighbors=&face_neighbors;
        if(face_index==0){helper.cell_to_add=0;MAP_BINTREE_MESH<T>::Map_Face_Process_Face(&helper,neighbor,this,Get_All_Face_Neighbors_Helper);}
        else{assert(face_index==1);helper.cell_to_add=1;MAP_BINTREE_MESH<T>::Map_Face_Process_Face(&helper,this,neighbor,Get_All_Face_Neighbors_Helper);}}
    else face_neighbors.Append(neighbor);
}
//#####################################################################
template class BINTREE_CELL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BINTREE_CELL<double>;
#endif
#endif
