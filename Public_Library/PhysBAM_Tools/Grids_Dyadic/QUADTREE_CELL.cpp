//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Cell_Compaction_Mapping
//#####################################################################
template<class T> void QUADTREE_CELL<T>::
Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count)
{
    mapping_array(Cell())=++cell_count;Cell()=cell_count;
    if(Has_Children()) for(int i=0;i<4;i++) Child(i)->Create_Cell_Compaction_Mapping(mapping_array,cell_count);
}
//#####################################################################
// Function Storage_Requirement
//#####################################################################
template<class T> int QUADTREE_CELL<T>::
Storage_Requirement() const
{
    int total=0;
    if(Has_Children()){total+=sizeof(QUADTREE_CHILDREN<T>);for(int i=0;i<4;i++) total+=Child(i)->Storage_Requirement();}
    return total;   
}
//#####################################################################
// Function Create_Root
//#####################################################################
template<class T> QUADTREE_CELL<T>* QUADTREE_CELL<T>::
Create_Root(const bool use_nodes,const bool use_faces,int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,2>& root_center,const VECTOR<T,2>& root_DX)
{
    QUADTREE_CHILDREN<T>* owner=new QUADTREE_CHILDREN<T>(use_nodes,use_faces);
    owner->Initialize_Root_Cell(number_of_cells,number_of_nodes,number_of_faces,root_center,root_DX);
    return &(owner->children[0]);
}
//#####################################################################
// Function Create_Children
//#####################################################################
template<class T> void QUADTREE_CELL<T>::
Create_Children(int& number_of_cells,ARRAY<QUADTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_nodes,int& number_of_faces,
                ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_faces,const QUADTREE_GRID<T>* grid)
{
    assert(children == 0);children=new QUADTREE_CHILDREN<T>(owner->nodes!=0,owner->faces!=0);
    children->Initialize_Non_Root_Cell(this,number_of_cells,new_cells,number_of_nodes,new_nodes,number_of_faces,new_faces,Center(),grid);
}
//#####################################################################
// Function Delete_Children
//#####################################################################
// returns an array of the nodes that were deleted
template<class T> void QUADTREE_CELL<T>::
Delete_Children()
{
    if(!Has_Children()) return;
    for(int i=0;i<4;i++) Child(i)->Delete_Children();
    delete children;children=0;
}
//#####################################################################
// Function Get_Neighbor
//#####################################################################
// returns the requested neighbor or 0 if no neighbor exists - offsets are between -1 and 1
// same_depth_only=true makes the function return a cell at the same depth or 0 - same_depth_only=false may return a neighbor cell higher up
template<class T> QUADTREE_CELL<T>* QUADTREE_CELL<T>::
Get_Neighbor(const int x_offset,const int y_offset,const QUADTREE_GRID<T>* grid,const ARRAY<VECTOR<QUADTREE_CELL<T>*,4> >* neighbors,const bool same_depth_only) const
{
    TV_INT offset(x_offset,y_offset);
    assert(RANGE<TV_INT>(-TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()).Lazy_Inside(offset));
    if(neighbors){
        assert(x_offset==0||y_offset==0);
        if(x_offset!=0) return (*neighbors)(Cell())(((x_offset+1)>>1)+1);
        assert(y_offset!=0);return (*neighbors)(Cell())(((y_offset+1)>>1)+3);}
    if(owner->parent == 0){ // we are at the top of the tree... use the uniform grid
        if(grid==0)return 0;TV_INT i=grid->uniform_grid.Index(Center());
        if(grid->cells.Domain_Indices().Lazy_Outside(i+offset)) return 0;
        return grid->cells(i+offset);}
    int possible_child=child_index+x_offset*X_SIBLING+y_offset*Y_SIBLING;
    if(Can_Go_In_X_Direction(child_index,x_offset) && Can_Go_In_Y_Direction(child_index,y_offset)){
        assert(possible_child >= 0 && possible_child < 4);return &owner->children[possible_child];}
    else{ // we need to keep going up the tree to find the node where we can go down the right branch
        if(x_offset!=0 && Can_Go_In_X_Direction(child_index,x_offset)){
            QUADTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(0,y_offset,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
            if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
            return correct_branch->Child(child_index+x_offset*X_SIBLING-y_offset*Y_SIBLING);}
        if(y_offset!=0 && Can_Go_In_Y_Direction(child_index, y_offset)){
            QUADTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(x_offset,0,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
            if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
            return correct_branch->Child(child_index-x_offset*X_SIBLING+y_offset*Y_SIBLING);}
        else{
            QUADTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(x_offset,y_offset,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
            if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
            return correct_branch->Child(child_index-x_offset*X_SIBLING-y_offset*Y_SIBLING);}}
}
//#####################################################################
// Function Can_Go_In_X_Direction
//#####################################################################
template<class T> bool QUADTREE_CELL<T>::
Can_Go_In_X_Direction(const int child_index,const int x_offset)
{
    if(x_offset == 0) return true;
    assert(x_offset >= -1 && x_offset <= 1);
    static const bool table[4]={true,false,true,false};
    if(x_offset == -1) return !table[child_index];else return table[child_index];
}
//#####################################################################
// Function Can_Go_In_Y_Direction
//#####################################################################
template<class T> bool QUADTREE_CELL<T>::
Can_Go_In_Y_Direction(const int child_index,const int y_offset)
{
    if(y_offset == 0)  return true;
    assert(y_offset >= -1 && y_offset <= 1);
    static const bool table[4]={true,true,false,false};
    if(y_offset == -1) return !table[child_index];else return table[child_index];
}
//#####################################################################
// Function Get_All_Face_Neighbors
//#####################################################################
template<class T> struct GET_ALL_FACE_NEIGHBORS_HELPER{int cell_to_add;ARRAY<QUADTREE_CELL<T>*>* face_neighbors;};
template<class T> static void
Get_All_Face_Neighbors_Helper(void* data,const QUADTREE_CELL<T>* cell1,const QUADTREE_CELL<T>* cell2,const int axis)
{
    GET_ALL_FACE_NEIGHBORS_HELPER<T>* helper=(GET_ALL_FACE_NEIGHBORS_HELPER<T>*)data;
    const QUADTREE_CELL<T>* cells[]={cell1,cell2};
    helper->face_neighbors->Append((QUADTREE_CELL<T>*)cells[helper->cell_to_add]);
}
//#####################################################################
// Function Get_All_Face_Neighbors
//#####################################################################
template<class T> void QUADTREE_CELL<T>::
Get_All_Face_Neighbors(int face_index,ARRAY<QUADTREE_CELL<T>*>& face_neighbors,const QUADTREE_GRID<T>* grid) const
{
    int face_index_to_offset[][2]={{-1,0},{1,0},{0,-1},{0,1}}; // convertes from face index to x and y offsets for use in neighbor search
    QUADTREE_CELL<T>* neighbor=Get_Neighbor(face_index_to_offset[face_index][0],face_index_to_offset[face_index][1],grid,0,false);
    if(!neighbor) return;
    if(neighbor->Has_Children()){
        GET_ALL_FACE_NEIGHBORS_HELPER<T> helper;helper.face_neighbors=&face_neighbors;
        if(face_index==0){helper.cell_to_add=0;MAP_QUADTREE_MESH<T>::Map_Face_Process_Face(&helper,neighbor,this,X_AXIS,Get_All_Face_Neighbors_Helper);}
        else if(face_index==1){helper.cell_to_add=1;MAP_QUADTREE_MESH<T>::Map_Face_Process_Face(&helper,this,neighbor,X_AXIS,Get_All_Face_Neighbors_Helper);}
        else if(face_index==2){helper.cell_to_add=0;MAP_QUADTREE_MESH<T>::Map_Face_Process_Face(&helper,neighbor,this,Y_AXIS,Get_All_Face_Neighbors_Helper);}
        else{assert(face_index==3);helper.cell_to_add=1;
            MAP_QUADTREE_MESH<T>::Map_Face_Process_Face(&helper,this,neighbor,Y_AXIS,Get_All_Face_Neighbors_Helper);}}
    else face_neighbors.Append(neighbor);
}
//#####################################################################
// Function Interpolate_Node_Values_To_Direct_Children
//#####################################################################
template<class T> template<class TV> void QUADTREE_CELL<T>::
Interpolate_Node_Values_To_Direct_Children(ARRAY<TV>& node_values)
{
    assert(Has_Children());
    node_values(children->Node(0,3))=(T).25*(node_values(Node(0))+node_values(Node(1))+node_values(Node(2))+node_values(Node(3)));
    node_values(children->Node(0,1))=(T).5*(node_values(Node(0))+node_values(Node(1)));
    node_values(children->Node(0,2))=(T).5*(node_values(Node(0))+node_values(Node(2)));
    node_values(children->Node(1,3))=(T).5*(node_values(Node(1))+node_values(Node(3)));
    node_values(children->Node(2,3))=(T).5*(node_values(Node(2))+node_values(Node(3)));
}
//#####################################################################
// Function Interpolate_Face_Values_To_Direct_Children
//#####################################################################
template<class T> template<class TV> void QUADTREE_CELL<T>::
Interpolate_Face_Values_To_Direct_Children(ARRAY<TV>& face_values)
{
    assert(Has_Children());
    face_values(children->Face(0,0))=face_values(children->Face(2,0))=face_values(Face(0));
    face_values(children->Face(1,1))=face_values(children->Face(3,1))=face_values(Face(1));
    face_values(children->Face(0,2))=face_values(children->Face(1,2))=face_values(Face(2));
    face_values(children->Face(2,3))=face_values(children->Face(3,3))=face_values(Face(3));
    face_values(children->Face(0,3))=face_values(children->Face(1,3))=(T).5*(face_values(Face(2))+face_values(Face(3)));
    face_values(children->Face(0,1))=face_values(children->Face(2,1))=(T).5*(face_values(Face(0))+face_values(Face(1)));
}
//#####################################################################
template class QUADTREE_CELL<float>;
template void QUADTREE_CELL<float>::Interpolate_Face_Values_To_Direct_Children<float>(ARRAY<float,int>&);
template void QUADTREE_CELL<float>::Interpolate_Node_Values_To_Direct_Children<float>(ARRAY<float,int>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class QUADTREE_CELL<double>;
template void QUADTREE_CELL<double>::Interpolate_Face_Values_To_Direct_Children<double>(ARRAY<double,int>&);
template void QUADTREE_CELL<double>::Interpolate_Node_Values_To_Direct_Children<double>(ARRAY<double,int>&);
#endif
#endif
