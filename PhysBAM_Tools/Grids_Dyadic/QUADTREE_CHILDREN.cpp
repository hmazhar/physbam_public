//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <cstring>
using namespace PhysBAM;
//#####################################################################
template<class T> const int QUADTREE_CHILDREN<T>::
node_indices[4][4] ={{0,1,3,4},{1,2,4,5},{3,4,6,7},{4,5,7,8}};
//#####################################################################
template<class T> const int QUADTREE_CHILDREN<T>::
face_indices[4][4] ={{0,1,6,7},{1,2,9,10},{3,4,7,8},{4,5,10,11}};
//#####################################################################
// Function Create_Face_Compaction_Mapping_Helper
//#####################################################################
template<class T> void QUADTREE_CHILDREN<T>::
Create_Face_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& face_count)
{
    for(int i=0;i<12;i++){
        if(mapping_array(faces[i])<=0) mapping_array(faces[i])=++face_count;
        faces[i]=mapping_array(faces[i]);}
    for(int i=0;i<4;i++) if(children[i].Has_Children()) children[i].children->Create_Face_Compaction_Mapping_Helper(mapping_array,face_count);
}
//#####################################################################
// Function Create_Node_Compaction_Mapping_Helper
//#####################################################################
template<class T> void QUADTREE_CHILDREN<T>::
Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count)
{
    for(int i=0;i<9;i++){
        if(mapping_array(nodes[i])<=0) mapping_array(nodes[i])=++node_count;
        nodes[i]=mapping_array(nodes[i]);}
    for(int i=0;i<4;i++) if(children[i].Has_Children()) children[i].children->Create_Node_Compaction_Mapping_Helper(mapping_array,node_count);
}
//#####################################################################
// Function Initialize_Root
//#####################################################################
template<class T> void QUADTREE_CHILDREN<T>::
Initialize_Root_Cell(int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,2>& root_center,const VECTOR<T,2>& root_DX)
{
    parent=0;childrens_depth=1;childrens_DX=root_DX;parents_center=root_center+(T).5*root_DX;
    int cell_input=0;children[0].Initialize(this,number_of_cells,cell_input);
    if(nodes){
        for(int i=0;i<9;i++) nodes[i]=-1;
        for(int i=0;i<4;i++) nodes[node_indices[0][i]]=++number_of_nodes;}
    if(faces){
        for(int i=0;i<12;i++) faces[i]=-1;
        for(int i=0;i<4;i++) faces[face_indices[0][i]]=++number_of_faces;}
}
//#####################################################################
// Function Initialize_Pseudo_Root_Cells
//#####################################################################
template<class T> void QUADTREE_CHILDREN<T>::
Initialize_Pseudo_Root_Cells(int& number_of_cells,ARRAY<int,VECTOR<int,1> >& nodes_input,ARRAY<int,VECTOR<int,1> >& faces_input,const VECTOR<T,2>& center_input,const VECTOR<T,2>& childrens_DX_input)
{
    parent=0;
    childrens_depth=1;
    parents_center=center_input;childrens_DX=childrens_DX_input;
    for(int i=0;i<4;i++) children[i].Initialize(this,number_of_cells,i);

    if(nodes)for(int i=0;i<9;i++)nodes[i]=nodes_input(i);
    if(faces)for(int i=0;i<12;i++)faces[i]=faces_input(i);
}
//#####################################################################
// Function Initialize_Non_Root_Cell
//#####################################################################
template<class T> void QUADTREE_CHILDREN<T>::
Initialize_Non_Root_Cell(QUADTREE_CELL<T>* parent_input,int& number_of_cells,ARRAY<QUADTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_nodes,
                         int& number_of_faces,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_faces,const VECTOR<T,2>& parents_center_in,const QUADTREE_GRID<T>* grid)
{
    parent=parent_input;
    childrens_depth=parent->owner->childrens_depth+1;
    parents_center=parents_center_in;
    childrens_DX=(T).5*parent->owner->childrens_DX;
    for(int i=0;i<4;i++){children[i].Initialize(this,number_of_cells,i);if(new_cells) new_cells->Append(&(children[i]));}

    QUADTREE_CELL<T> *cell_left,*cell_right,*cell_bottom,*cell_top;
    QUADTREE_CHILDREN<T> *children_left=0,*children_right=0,*children_bottom=0,*children_top=0;
    
    if(nodes || faces){ // get all the cells that are needed
        cell_left=parent->Get_Neighbor(-1,0,grid);cell_right=parent->Get_Neighbor(1,0,grid);
        cell_bottom=parent->Get_Neighbor(0,-1,grid);cell_top=parent->Get_Neighbor(0,1,grid);
        children_left=(cell_left==0)?0:cell_left->children;children_right=(cell_right==0)?0:cell_right->children;
        children_bottom=(cell_bottom==0)?0:cell_bottom->children;children_top=(cell_top==0)?0:cell_top->children;}
    
    if(nodes){
        memset(nodes,0,sizeof(int)*9);
        // corner nodes from the parent cell
        nodes[Node_Index(0,0)]=parent->Node(0);nodes[Node_Index(2,0)]=parent->Node(1);nodes[Node_Index(0,2)]=parent->Node(2);nodes[Node_Index(2,2)]=parent->Node(3);
        // add a new center node
        nodes[Node_Index(1,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[0],3));
        if(children_left) nodes[Node_Index(0,1)]=children_left->nodes[Node_Index(2,1)];
        else{nodes[Node_Index(0,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[0],2));}
        if(children_right) nodes[Node_Index(2,1)]=children_right->nodes[Node_Index(0,1)];
        else{nodes[Node_Index(2,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[1],3));}
        if(children_bottom) nodes[Node_Index(1,0)]=children_bottom->nodes[Node_Index(1,2)];
        else{nodes[Node_Index(1,0)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[0],1));}
        if(children_top) nodes[Node_Index(1,2)]=children_top->nodes[Node_Index(1,0)];
        else{nodes[Node_Index(1,2)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[2],3));}}

    if(faces){ // add the new internal faces
        faces[face_indices[0][1]]=++number_of_faces;faces[face_indices[2][1]]=++number_of_faces; // x-faces
        faces[face_indices[0][3]]=++number_of_faces;faces[face_indices[1][3]]=++number_of_faces; // y-faces
        if(new_faces){
            new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[0],1));new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[2],1));
            new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[0],3));new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[1],3));}
        if(children_left){
            faces[face_indices[0][0]]=children_left->faces[face_indices[1][1]];faces[face_indices[2][0]]=children_left->faces[face_indices[3][1]];}
        else{
            faces[face_indices[0][0]]=++number_of_faces;faces[face_indices[2][0]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[0],0));new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[2],0));}}
        if(children_right){
            faces[face_indices[1][1]]=children_right->faces[face_indices[0][0]];faces[face_indices[3][1]]=children_right->faces[face_indices[2][0]];}
        else{
            faces[face_indices[1][1]]=++number_of_faces;faces[face_indices[3][1]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[1],1));new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[3],1));}}
        if(children_bottom){
            faces[face_indices[0][2]]=children_bottom->faces[face_indices[2][3]];faces[face_indices[1][2]]=children_bottom->faces[face_indices[3][3]];}
        else{
            faces[face_indices[0][2]]=++number_of_faces;faces[face_indices[1][2]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[0],2));new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[1],2));}}
        if(children_top){
            faces[face_indices[2][3]]=children_top->faces[face_indices[0][2]];faces[face_indices[3][3]]=children_top->faces[face_indices[1][2]];}
        else{
            faces[face_indices[2][3]]=++number_of_faces;faces[face_indices[3][3]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[2],3));new_faces->Append(PAIR<QUADTREE_CELL<T>*,int>(&children[3],3));}}}
}
//#####################################################################
template class QUADTREE_CHILDREN<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class QUADTREE_CHILDREN<double>;
#endif
#endif
