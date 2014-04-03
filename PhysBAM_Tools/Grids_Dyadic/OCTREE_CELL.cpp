//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OCTREE_CELL
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Cell_Compaction_Mapping
//#####################################################################
template<class T> void OCTREE_CELL<T>::
Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count)
{
    mapping_array(Cell())=++cell_count;Cell()=cell_count;
    if(Has_Children()) for(int i=0;i<8;i++) Child(i)->Create_Cell_Compaction_Mapping(mapping_array,cell_count);
}
//#####################################################################
// Function Storage_Requirement
//#####################################################################
template<class T> int OCTREE_CELL<T>::
Storage_Requirement() const
{
    int total=sizeof(OCTREE_CELL<T>);
    if(Has_Children()) total+=children->Storage_Requirement();
    return total;   
}
//#####################################################################
// Function Create_Root
//#####################################################################
template<class T> OCTREE_CELL<T>* OCTREE_CELL<T>::
Create_Root(const bool use_nodes,const bool use_faces,int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,3>& root_center,const VECTOR<T,3>& root_DX)
{
    OCTREE_CHILDREN<T>* owner=new OCTREE_CHILDREN<T>(use_nodes,use_faces);
    owner->Initialize_Root_Cell(number_of_cells,number_of_nodes,number_of_faces,root_center,root_DX);
    return &(owner->children[0]);
}
//#####################################################################
// Function Create_Children
//#####################################################################
template<class T> void OCTREE_CELL<T>::
Create_Children(int& number_of_cells,ARRAY<OCTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes,int& number_of_faces,
                ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_faces,const OCTREE_GRID<T>* grid)
{
    assert(children == 0);children=new OCTREE_CHILDREN<T>(owner->nodes!=0,owner->faces!=0);
    children->Initialize_Non_Root_Cell(this,number_of_cells,new_cells,number_of_nodes,new_nodes,number_of_faces,new_faces,Center(),grid);
}
//#####################################################################
// Function Delete_Children
//#####################################################################
// returns an array of the nodes that were deleted
template<class T> void OCTREE_CELL<T>::
Delete_Children()
{
    if(!Has_Children()) return;
    for(int i=0;i<8;i++) Child(i)->Delete_Children();
    delete children;children=0;
}
//#####################################################################
// Function Get_Neighbor
//#####################################################################
// returns the requested neighbor or 0 if no neighbor exists - offsets are between -1 and 1
// same_depth_only=true makes the function return a cell at the same depth or 0 - same_depth_only=false may return a neighbor cell higher up
template<class T> OCTREE_CELL<T>* OCTREE_CELL<T>::
Get_Neighbor(const int x_offset,const int y_offset,const int z_offset,const OCTREE_GRID<T>* grid,const ARRAY<VECTOR<OCTREE_CELL<T>*,6> >* neighbors,const bool same_depth_only) const
{
    TV_INT offset(x_offset,y_offset,z_offset);
    assert(RANGE<TV_INT>(-TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()).Lazy_Inside(offset));
    assert(offset.Contains(0));
    if(neighbors){
        assert((x_offset==0&&y_offset==0) || (x_offset==0&&z_offset==0) || (y_offset==0&&z_offset==0));
        if(x_offset!=0) return (*neighbors)(cell)(((x_offset+1)>>1)+1);
        if(y_offset!=0) return (*neighbors)(cell)(((y_offset+1)>>1)+3);
        assert(z_offset!=0);return (*neighbors)(cell)(((z_offset+1)>>1)+5);}
    if(owner->parent == 0){ // we are at the top of the tree... use the uniform uniform_grid
        if(grid==0)return 0;TV_INT i=grid->uniform_grid.Index(Center());
        if(grid->cells.Domain_Indices().Lazy_Outside(i+offset)) return 0;
        return grid->cells(i+offset);}
    int possible_child=child_index+x_offset*X_SIBLING+y_offset*Y_SIBLING+z_offset*Z_SIBLING;
    if(Can_Go_In_X_Direction(child_index,x_offset) && Can_Go_In_Y_Direction(child_index,y_offset) && Can_Go_In_Z_Direction(child_index,z_offset)){
        assert(possible_child >= 0 && possible_child < 8);return &owner->children[possible_child];}
    else{ // we need to keep going up the tree to find the node where we can go down the right branch
        if(x_offset!=0 && Can_Go_In_X_Direction(child_index,x_offset)){
            OCTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(0,y_offset,z_offset,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
            if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
            return correct_branch->Child(child_index+x_offset*X_SIBLING-y_offset*Y_SIBLING-z_offset*Z_SIBLING);}
        if(y_offset!=0 && Can_Go_In_Y_Direction(child_index, y_offset)){
            OCTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(x_offset,0,z_offset,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
            if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
            return correct_branch->Child(child_index-x_offset*X_SIBLING+y_offset*Y_SIBLING-z_offset*Z_SIBLING);}
        if(z_offset!=0 && Can_Go_In_Z_Direction(child_index,z_offset)){
            OCTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(x_offset,y_offset,0,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
            if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
            return correct_branch->Child(child_index-x_offset*X_SIBLING-y_offset*Y_SIBLING+z_offset*Z_SIBLING);}
        else{
            OCTREE_CELL<T>* correct_branch=owner->parent->Get_Neighbor(x_offset,y_offset,z_offset,grid,neighbors,same_depth_only);if(correct_branch == 0) return 0;
            if(!correct_branch->Has_Children()){if(same_depth_only) return 0;else return correct_branch;}
            return correct_branch->Child(child_index-x_offset*X_SIBLING-y_offset*Y_SIBLING-z_offset*Z_SIBLING);}}
}
//#####################################################################
// Function Can_Go_In_X_Direction
//#####################################################################
template<class T> bool OCTREE_CELL<T>::
Can_Go_In_X_Direction(const int child_index,const int x_offset)
{
    if(x_offset == 0) return true;
    assert(x_offset >= -1 && x_offset <= 1);
    static const bool table[8]={true,false,true,false,true,false,true,false};
    if(x_offset == -1) return !table[child_index];else return table[child_index];
}
//#####################################################################
// Function Can_Go_In_Y_Direction
//#####################################################################
template<class T> bool OCTREE_CELL<T>::
Can_Go_In_Y_Direction(const int child_index,const int y_offset)
{
    if(y_offset == 0)  return true;
    assert(y_offset >= -1 && y_offset <= 1);
    static const bool table[8]={true,true,false,false,true,true,false,false};
    if(y_offset == -1) return !table[child_index];else return table[child_index];
}
//#####################################################################
// Function Can_Go_In_Z_Direction
//#####################################################################
template<class T> bool OCTREE_CELL<T>::
Can_Go_In_Z_Direction(const int child_index,const int z_offset)
{
    if(z_offset == 0)  return true;
    assert(z_offset >= -1 && z_offset <= 1);
    static const bool table[8]={true,true,true,true,false,false,false,false};
    if(z_offset == -1) return !table[child_index];else return table[child_index];
}
//#####################################################################
// Function Get_All_Face_Neighbors
//#####################################################################
namespace PhysBAM{template<class T> struct GET_ALL_FACE_NEIGHBORS_HELPER{int cell_to_add;ARRAY<OCTREE_CELL<T>*>* face_neighbors;};}
template<class T> static void
Get_All_Face_Neighbors_Helper(void* data,const OCTREE_CELL<T>* cell1,const OCTREE_CELL<T>* cell2,const int axis)
{
    GET_ALL_FACE_NEIGHBORS_HELPER<T>* helper=(GET_ALL_FACE_NEIGHBORS_HELPER<T>*)data;
    const OCTREE_CELL<T>* cells[]={cell1,cell2};
    helper->face_neighbors->Append((OCTREE_CELL<T>*)cells[helper->cell_to_add]);
}
template<class T> void OCTREE_CELL<T>::
Get_All_Face_Neighbors(int face_index,ARRAY<OCTREE_CELL<T>*>& face_neighbors,const OCTREE_GRID<T>* grid) const
{
    int face_index_to_offset[][3]={{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}}; // convertes from face index to x and y offsets for use in neighbor search
    OCTREE_CELL<T>* neighbor=Get_Neighbor(face_index_to_offset[face_index][0],face_index_to_offset[face_index][1],face_index_to_offset[face_index][2],grid,0,false);
    if(!neighbor) return;
    if(neighbor->Has_Children()){
        GET_ALL_FACE_NEIGHBORS_HELPER<T> helper;helper.face_neighbors=&face_neighbors;
        if(face_index==0){helper.cell_to_add=0;MAP_OCTREE_MESH<T>::Map_Face_Process_Face(&helper,neighbor,this,X_AXIS,Get_All_Face_Neighbors_Helper);}
        else if(face_index==1){helper.cell_to_add=1;MAP_OCTREE_MESH<T>::Map_Face_Process_Face(&helper,this,neighbor,X_AXIS,Get_All_Face_Neighbors_Helper);}
        else if(face_index==2){helper.cell_to_add=0;MAP_OCTREE_MESH<T>::Map_Face_Process_Face(&helper,neighbor,this,Y_AXIS,Get_All_Face_Neighbors_Helper);}
        else if(face_index==3){helper.cell_to_add=1;MAP_OCTREE_MESH<T>::Map_Face_Process_Face(&helper,this,neighbor,Y_AXIS,Get_All_Face_Neighbors_Helper);}
        else if(face_index==4){helper.cell_to_add=0;MAP_OCTREE_MESH<T>::Map_Face_Process_Face(&helper,neighbor,this,Z_AXIS,Get_All_Face_Neighbors_Helper);}
        else{assert(face_index==5);helper.cell_to_add=1;MAP_OCTREE_MESH<T>::Map_Face_Process_Face(&helper,this,neighbor,Z_AXIS,Get_All_Face_Neighbors_Helper);}}
    else face_neighbors.Append(neighbor);
}
//#####################################################################
// Function Interpolate_Node_Values_To_Direct_Children
//#####################################################################
template<class T> template<class TV> void OCTREE_CELL<T>::
Interpolate_Node_Values_To_Direct_Children(ARRAY<TV>& node_values)
{
    assert(Has_Children());
    node_values(children->Node(0,7))=(T).125*(node_values(Node(0))+node_values(Node(1))+node_values(Node(2))+node_values(Node(3))+node_values(Node(4))+node_values(Node(5))+node_values(Node(6))
        +node_values(Node(7))); // middle
    node_values(children->Node(0,6))=(T).25*(node_values(Node(0))+node_values(Node(2))+node_values(Node(4))+node_values(Node(6))); // left
    node_values(children->Node(1,7))=(T).25*(node_values(Node(1))+node_values(Node(3))+node_values(Node(5))+node_values(Node(7))); // right
    node_values(children->Node(0,5))=(T).25*(node_values(Node(0))+node_values(Node(1))+node_values(Node(4))+node_values(Node(5))); // bottom
    node_values(children->Node(2,7))=(T).25*(node_values(Node(2))+node_values(Node(3))+node_values(Node(6))+node_values(Node(7))); // top
    node_values(children->Node(0,3))=(T).25*(node_values(Node(0))+node_values(Node(1))+node_values(Node(2))+node_values(Node(3))); // front
    node_values(children->Node(4,7))=(T).25*(node_values(Node(4))+node_values(Node(5))+node_values(Node(6))+node_values(Node(7))); // back
    node_values(children->Node(0,4))=(T).5*(node_values(Node(0))+node_values(Node(4)));node_values(children->Node(1,5))=(T).5*(node_values(Node(1))+node_values(Node(5)));
    node_values(children->Node(2,6))=(T).5*(node_values(Node(2))+node_values(Node(6)));node_values(children->Node(3,7))=(T).5*(node_values(Node(3))+node_values(Node(7)));
    node_values(children->Node(0,1))=(T).5*(node_values(Node(0))+node_values(Node(1)));node_values(children->Node(4,5))=(T).5*(node_values(Node(4))+node_values(Node(5)));
    node_values(children->Node(2,3))=(T).5*(node_values(Node(2))+node_values(Node(3)));node_values(children->Node(6,7))=(T).5*(node_values(Node(6))+node_values(Node(7)));
    node_values(children->Node(0,2))=(T).5*(node_values(Node(0))+node_values(Node(2)));node_values(children->Node(1,3))=(T).5*(node_values(Node(1))+node_values(Node(3)));
    node_values(children->Node(4,6))=(T).5*(node_values(Node(4))+node_values(Node(6)));node_values(children->Node(5,7))=(T).5*(node_values(Node(5))+node_values(Node(7)));
}
//#####################################################################
// Function Interpolate_Face_Values_To_Direct_Children
//#####################################################################
template<class T> template<class TV> void OCTREE_CELL<T>::
Interpolate_Face_Values_To_Direct_Children(ARRAY<TV>& face_values)
{
    assert(Has_Children());
    face_values(children->Face(0,0))=face_values(children->Face(2,0))=face_values(children->Face(4,0))=face_values(children->Face(6,0))=face_values(Face(0)); // LX external
    face_values(children->Face(1,1))=face_values(children->Face(3,1))=face_values(children->Face(5,1))=face_values(children->Face(7,1))=face_values(Face(1)); // HX external
    face_values(children->Face(0,2))=face_values(children->Face(1,2))=face_values(children->Face(4,2))=face_values(children->Face(5,2))=face_values(Face(2)); // LY external
    face_values(children->Face(2,3))=face_values(children->Face(3,3))=face_values(children->Face(6,3))=face_values(children->Face(7,3))=face_values(Face(3)); // HY external
    face_values(children->Face(0,4))=face_values(children->Face(1,4))=face_values(children->Face(2,4))=face_values(children->Face(3,4))=face_values(Face(4)); // LZ external
    face_values(children->Face(4,5))=face_values(children->Face(5,5))=face_values(children->Face(6,5))=face_values(children->Face(7,5))=face_values(Face(5)); // HZ external
    face_values(children->Face(0,1))=face_values(children->Face(2,1))=face_values(children->Face(4,1))=face_values(children->Face(6,1))=(T).5*(face_values(Face(0))+face_values(Face(1))); // x internal
    face_values(children->Face(0,3))=face_values(children->Face(1,3))=face_values(children->Face(4,3))=face_values(children->Face(5,3))=(T).5*(face_values(Face(2))+face_values(Face(3))); // y internal
    face_values(children->Face(0,5))=face_values(children->Face(1,5))=face_values(children->Face(2,5))=face_values(children->Face(3,5))=(T).5*(face_values(Face(4))+face_values(Face(5))); // z internal
}
//#####################################################################
template class OCTREE_CELL<float>;
template void OCTREE_CELL<float>::Interpolate_Node_Values_To_Direct_Children<float>(ARRAY<float,int>&);
template void OCTREE_CELL<float>::Interpolate_Face_Values_To_Direct_Children<float>(ARRAY<float,int>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OCTREE_CELL<double>;
template void OCTREE_CELL<double>::Interpolate_Node_Values_To_Direct_Children<double>(ARRAY<double,int>&);
template void OCTREE_CELL<double>::Interpolate_Face_Values_To_Direct_Children<double>(ARRAY<double,int>&);
#endif
#endif
