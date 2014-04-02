//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <cstring>
using namespace PhysBAM;
//#####################################################################
template<class T> const int BINTREE_CHILDREN<T>::
node_indices[2][2] = {{0,1},{1,2}};
//#####################################################################
template<class T> const int BINTREE_CHILDREN<T>::
face_indices[2][2] = {{0,1},{1,2}};
//#####################################################################
// Function Create_Face_Compaction_Mapping_Helper
//#####################################################################
template<class T> void BINTREE_CHILDREN<T>::
Create_Face_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& face_count)
{
    for(int i=0;i<3;i++){
        if(mapping_array(faces[i])<=0) mapping_array(faces[i])=++face_count;
        faces[i]=mapping_array(faces[i]);}
    for(int i=0;i<2;i++) if(children[i].Has_Children()) children[i].children->Create_Face_Compaction_Mapping_Helper(mapping_array,face_count);
}
//#####################################################################
// Function Create_Node_Compaction_Mapping_Helper
//#####################################################################
template<class T> void BINTREE_CHILDREN<T>::
Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count)
{
    for(int i=0;i<3;i++){
        if(mapping_array(nodes[i])<=0) mapping_array(nodes[i])=++node_count;
        nodes[i]=mapping_array(nodes[i]);}
    for(int i=0;i<2;i++) if(children[i].Has_Children()) children[i].children->Create_Node_Compaction_Mapping_Helper(mapping_array,node_count);
}
//#####################################################################
// Function Initialize_Root
//#####################################################################
template<class T> void BINTREE_CHILDREN<T>::
Initialize_Root_Cell(int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,1>& root_center,const VECTOR<T,1>& root_DX)
{
    parent=0;childrens_depth=1;childrens_DX=root_DX;parents_center=root_center+(T).5*root_DX;
    int cell_input=0;children[0].Initialize(this,number_of_cells,cell_input);
    if(nodes){
        for(int i=0;i<3;i++) nodes[i]=-1;
        for(int i=0;i<2;i++) nodes[node_indices[0][i]]=++number_of_nodes;}
    if(faces){
        for(int i=0;i<3;i++) faces[i]=-1;
        for(int i=0;i<2;i++) faces[face_indices[0][i]]=++number_of_faces;}
}
//#####################################################################
// Function Initialize_Pseudo_Root_Cells
//#####################################################################
template<class T> void BINTREE_CHILDREN<T>::
Initialize_Pseudo_Root_Cells(int& number_of_cells,ARRAY<int,VECTOR<int,1> >& nodes_input,ARRAY<int,VECTOR<int,1> >& faces_input,const VECTOR<T,1>& center_input,const VECTOR<T,1>& childrens_DX_input)
{
    parent=0;
    childrens_depth=1;
    parents_center=center_input;childrens_DX=childrens_DX_input;
    for(int i=0;i<2;i++) children[i].Initialize(this,number_of_cells,i);

    if(nodes) for(int i=0;i<3;i++)nodes[i]=nodes_input(i);
    if(faces) for(int i=0;i<3;i++)faces[i]=faces_input(i);
}
//#####################################################################
// Function Initialize_Non_Root_Cell
//#####################################################################
template<class T> void BINTREE_CHILDREN<T>::
Initialize_Non_Root_Cell(BINTREE_CELL<T>* parent_input,int& number_of_cells,ARRAY<BINTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_nodes,
                         int& number_of_faces,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_faces,const VECTOR<T,1>& parents_center_in,const BINTREE_GRID<T>* grid)
{
    parent=parent_input;
    childrens_depth=parent->owner->childrens_depth+1;
    parents_center=parents_center_in;
    childrens_DX=(T).5*parent->owner->childrens_DX;
    for(int i=0;i<2;i++){children[i].Initialize(this,number_of_cells,i);if(new_cells) new_cells->Append(&(children[i]));}

    BINTREE_CELL<T> *cell_left,*cell_right;
    BINTREE_CHILDREN<T> *children_left=0,*children_right=0;
    
    if(nodes || faces){ // get all the cells that are needed
        cell_left=parent->Get_Neighbor(-1,grid);cell_right=parent->Get_Neighbor(1,grid);
        children_left=(cell_left==0)?0:cell_left->children;children_right=(cell_right==0)?0:cell_right->children;}
    
    if(nodes){
        memset(nodes,0,sizeof(int)*3);
        nodes[0]=parent->Node(0);nodes[2]=parent->Node(1); // corner nodes from the parent cell
        nodes[1]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<BINTREE_CELL<T>*,int>(&children[0],1));}

    if(faces){ // add the new internal faces
        faces[face_indices[0][0]]=parent->Face(0);
        faces[face_indices[1][1]]=parent->Face(1);
        faces[face_indices[0][1]]=++number_of_faces;
        if(new_faces) new_faces->Append(PAIR<BINTREE_CELL<T>*,int>(&children[0],1));}
}
//#####################################################################
template class BINTREE_CHILDREN<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BINTREE_CHILDREN<double>;
#endif
#endif
