//#####################################################################
// Copyright 2003-2004, Ron Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <cstring>
using namespace PhysBAM;
//#####################################################################
template<class T> const int OCTREE_CHILDREN<T>::
node_indices[8][8] ={{0,1,3,4,9,10,12,13},{1,2,4,5,10,11,13,14},{3,4,6,7,12,13,15,16},{4,5,7,8,13,14,16,17},{9,10,12,13,18,19,21,22},{10,11,13,14,19,20,22,23},{12,13,15,16,21,22,24,25},{13,14,16,17,22,23,25,26}};
//#####################################################################
template<class T> const int OCTREE_CHILDREN<T>::
face_indices[8][6] ={{0,1,12,13,24,25},{1,2,15,16,27,28},{3,4,13,14,30,31},{4,5,16,17,33,34},{6,7,18,19,25,26},{7,8,21,22,28,29},{9,10,19,20,31,32},{10,11,22,23,34,35}};
//#####################################################################
// Function Create_Face_Compaction_Mapping_Helper
//#####################################################################
template<class T> void OCTREE_CHILDREN<T>::
Create_Face_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& face_count)
{
    for(int i=0;i<36;i++){
        if(mapping_array(faces[i])<=0) mapping_array(faces[i])=++face_count;
        faces[i]=mapping_array(faces[i]);}
    for(int i=0;i<8;i++) if(children[i].Has_Children()) children[i].children->Create_Face_Compaction_Mapping_Helper(mapping_array,face_count);
}
//#####################################################################
// Function Create_Node_Compaction_Mapping_Helper
//#####################################################################
template<class T> void OCTREE_CHILDREN<T>::
Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count)
{
    for(int i=0;i<27;i++){
        if(mapping_array(nodes[i])<=0) mapping_array(nodes[i])=++node_count;
        nodes[i]=mapping_array(nodes[i]);}
    for(int i=0;i<8;i++) if(children[i].Has_Children()) children[i].children->Create_Node_Compaction_Mapping_Helper(mapping_array,node_count);
}
//#####################################################################
// Function Initialize_Root
//#####################################################################
template<class T> void OCTREE_CHILDREN<T>::
Initialize_Root_Cell(int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,3>& root_center,const VECTOR<T,3>& root_DX)
{
    parent=0;childrens_depth=1;childrens_DX=root_DX;parents_center=root_center+(T).5*root_DX;
    children[0].Initialize(this,number_of_cells,0);
    if(nodes){
        for(int i=0;i<27;i++) nodes[i]=-1;
        for(int i=0;i<8;i++) nodes[node_indices[0][i]]=++number_of_nodes;}
    if(faces){
        for(int i=0;i<36;i++) faces[i]=-1;
        for(int i=0;i<6;i++) faces[face_indices[0][i]]=++number_of_faces;}
}
//#####################################################################
// Function Initialize_Pseudo_Root_Cells
//#####################################################################
template<class T> void OCTREE_CHILDREN<T>::
Initialize_Pseudo_Root_Cells(int& number_of_cells,ARRAY<int,VECTOR<int,1> >& nodes_input,ARRAY<int,VECTOR<int,1> >& faces_input,const VECTOR<T,3>& center_input,const VECTOR<T,3>& childrens_DX_input)
{
    parent=0;
    childrens_depth=1;
    parents_center=center_input;childrens_DX=childrens_DX_input;
    for(int i=0;i<8;i++) children[i].Initialize(this,number_of_cells,i);

    if(nodes)for(int i=0;i<27;i++)nodes[i]=nodes_input(i);
    if(faces)for(int i=0;i<36;i++)faces[i]=faces_input(i);
}
//#####################################################################
// Function Initialize_Non_Root_Cell
//#####################################################################
template<class T> void OCTREE_CHILDREN<T>::
Initialize_Non_Root_Cell(OCTREE_CELL<T>* parent_input,int& number_of_cells,ARRAY<OCTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes,
                         int& number_of_faces,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_faces,const VECTOR<T,3>& parents_center_in,const OCTREE_GRID<T>* grid)
{
    parent=parent_input;
    childrens_depth=parent->owner->childrens_depth+1;
    parents_center=parents_center_in;
    childrens_DX=(T).5*parent->owner->childrens_DX;
    for(int i=0;i<8;i++){children[i].Initialize(this,number_of_cells,i);if(new_cells) new_cells->Append(&(children[i]));}

    OCTREE_CELL<T> *cell_left,*cell_right,*cell_bottom,*cell_top,*cell_front,*cell_back;
    OCTREE_CHILDREN<T> *children_left=0,*children_right=0,*children_bottom=0,*children_top=0,*children_front=0,*children_back=0;
    
    if(nodes || faces){ // get all the cells that are needed
        cell_left=parent->Get_Neighbor(-1,0,0,grid);cell_right=parent->Get_Neighbor(1,0,0,grid);
        cell_bottom=parent->Get_Neighbor(0,-1,0,grid);cell_top=parent->Get_Neighbor(0,1,0,grid);
        cell_front=parent->Get_Neighbor(0,0,-1,grid);cell_back=parent->Get_Neighbor(0,0,1,grid);
        children_left=(cell_left==0)?0:cell_left->children;children_right=(cell_right==0)?0:cell_right->children;
        children_bottom=(cell_bottom==0)?0:cell_bottom->children;children_top=(cell_top==0)?0:cell_top->children;
        children_front=(cell_front==0)?0:cell_front->children;children_back=(cell_back==0)?0:cell_back->children;}
    
    if(nodes){
        memset(nodes,0,sizeof(int)*27);
        // corner nodes from the parent cell
        nodes[Node_Index(0,0,0)]=parent->Node(0);nodes[Node_Index(2,0,0)]=parent->Node(1);nodes[Node_Index(0,2,0)]=parent->Node(2);nodes[Node_Index(2,2,0)]=parent->Node(3);
        nodes[Node_Index(0,0,2)]=parent->Node(4);nodes[Node_Index(2,0,2)]=parent->Node(5);nodes[Node_Index(0,2,2)]=parent->Node(6);nodes[Node_Index(2,2,2)]=parent->Node(7);
        // add a new center node
        nodes[Node_Index(1,1,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],7));
        if(children_left){
            nodes[Node_Index(0,0,1)]=children_left->nodes[Node_Index(2,0,1)];nodes[Node_Index(0,1,0)]=children_left->nodes[Node_Index(2,1,0)];
            nodes[Node_Index(0,1,1)]=children_left->nodes[Node_Index(2,1,1)];nodes[Node_Index(0,1,2)]=children_left->nodes[Node_Index(2,1,2)];
            nodes[Node_Index(0,2,1)]=children_left->nodes[Node_Index(2,2,1)];}
        else{nodes[Node_Index(0,1,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],6));}
        if(children_right){
            nodes[Node_Index(2,0,1)]=children_right->nodes[Node_Index(0,0,1)];nodes[Node_Index(2,1,0)]=children_right->nodes[Node_Index(0,1,0)];
            nodes[Node_Index(2,1,1)]=children_right->nodes[Node_Index(0,1,1)];nodes[Node_Index(2,1,2)]=children_right->nodes[Node_Index(0,1,2)];
            nodes[Node_Index(2,2,1)]=children_right->nodes[Node_Index(0,2,1)];}
        else{nodes[Node_Index(2,1,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[1],7));}
        if(children_bottom){
            nodes[Node_Index(0,0,1)]=children_bottom->nodes[Node_Index(0,2,1)];nodes[Node_Index(1,0,0)]=children_bottom->nodes[Node_Index(1,2,0)];
            nodes[Node_Index(1,0,1)]=children_bottom->nodes[Node_Index(1,2,1)];nodes[Node_Index(1,0,2)]=children_bottom->nodes[Node_Index(1,2,2)];
            nodes[Node_Index(2,0,1)]=children_bottom->nodes[Node_Index(2,2,1)];}
        else{nodes[Node_Index(1,0,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],5));}
        if(children_top){
            nodes[Node_Index(0,2,1)]=children_top->nodes[Node_Index(0,0,1)];nodes[Node_Index(1,2,0)]=children_top->nodes[Node_Index(1,0,0)];
            nodes[Node_Index(1,2,1)]=children_top->nodes[Node_Index(1,0,1)];nodes[Node_Index(1,2,2)]=children_top->nodes[Node_Index(1,0,2)];
            nodes[Node_Index(2,2,1)]=children_top->nodes[Node_Index(2,0,1)];}
        else{nodes[Node_Index(1,2,1)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[2],7));}
        if(children_front){
            nodes[Node_Index(0,1,0)]=children_front->nodes[Node_Index(0,1,2)];nodes[Node_Index(1,0,0)]=children_front->nodes[Node_Index(1,0,2)];
            nodes[Node_Index(1,1,0)]=children_front->nodes[Node_Index(1,1,2)];nodes[Node_Index(1,2,0)]=children_front->nodes[Node_Index(1,2,2)];
            nodes[Node_Index(2,1,0)]=children_front->nodes[Node_Index(2,1,2)];}
        else{nodes[Node_Index(1,1,0)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],3));}
        if(children_back){
            nodes[Node_Index(0,1,2)]=children_back->nodes[Node_Index(0,1,0)];nodes[Node_Index(1,0,2)]=children_back->nodes[Node_Index(1,0,0)];
            nodes[Node_Index(1,1,2)]=children_back->nodes[Node_Index(1,1,0)];nodes[Node_Index(1,2,2)]=children_back->nodes[Node_Index(1,2,0)];
            nodes[Node_Index(2,1,2)]=children_back->nodes[Node_Index(2,1,0)];}
        else{nodes[Node_Index(1,1,2)]=++number_of_nodes;if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[4],7));}
        // check twelve edge neighbors for existing nodes
        Edge_Neighbor_Check(number_of_nodes,-1,-1,0,Node_Index(0,0,1),Node_Index(2,2,1),0,4,new_nodes,grid);Edge_Neighbor_Check(number_of_nodes,1,-1,0,Node_Index(2,0,1),Node_Index(0,2,1),1,5,new_nodes,grid);
        Edge_Neighbor_Check(number_of_nodes,-1,1,0,Node_Index(0,2,1),Node_Index(2,0,1),2,6,new_nodes,grid);Edge_Neighbor_Check(number_of_nodes,1,1,0,Node_Index(2,2,1),Node_Index(0,0,1),3,7,new_nodes,grid);
        Edge_Neighbor_Check(number_of_nodes,-1,0,-1,Node_Index(0,1,0),Node_Index(2,1,2),0,2,new_nodes,grid);Edge_Neighbor_Check(number_of_nodes,1,0,-1,Node_Index(2,1,0),Node_Index(0,1,2),1,3,new_nodes,grid);
        Edge_Neighbor_Check(number_of_nodes,-1,0,1,Node_Index(0,1,2),Node_Index(2,1,0),4,6,new_nodes,grid);Edge_Neighbor_Check(number_of_nodes,1,0,1,Node_Index(2,1,2),Node_Index(0,1,0),5,7,new_nodes,grid);
        Edge_Neighbor_Check(number_of_nodes,0,-1,-1,Node_Index(1,0,0),Node_Index(1,2,2),0,1,new_nodes,grid);Edge_Neighbor_Check(number_of_nodes,0,1,-1,Node_Index(1,2,0),Node_Index(1,0,2),2,3,new_nodes,grid);
        Edge_Neighbor_Check(number_of_nodes,0,-1,1,Node_Index(1,0,2),Node_Index(1,2,0),4,5,new_nodes,grid);Edge_Neighbor_Check(number_of_nodes,0,1,1,Node_Index(1,2,2),Node_Index(1,0,0),6,7,new_nodes,grid);}

    if(faces){ // add the new internal faces
        faces[face_indices[0][1]]=++number_of_faces;faces[face_indices[2][1]]=++number_of_faces;faces[face_indices[4][1]]=++number_of_faces;faces[face_indices[6][1]]=++number_of_faces; // x-faces
        faces[face_indices[0][3]]=++number_of_faces;faces[face_indices[1][3]]=++number_of_faces;faces[face_indices[4][3]]=++number_of_faces;faces[face_indices[5][3]]=++number_of_faces; // y-faces
        faces[face_indices[0][5]]=++number_of_faces;faces[face_indices[1][5]]=++number_of_faces;faces[face_indices[2][5]]=++number_of_faces;faces[face_indices[3][5]]=++number_of_faces; // z-faces
        if(new_faces){
            new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],1));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[2],1));
            new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[4],1));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[6],1));
            new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],3));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[1],3));
            new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[4],3));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[5],3));
            new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],5));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[1],5));
            new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[2],5));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[3],5));}
        if(children_left){
            faces[face_indices[0][0]]=children_left->faces[face_indices[1][1]];faces[face_indices[2][0]]=children_left->faces[face_indices[3][1]];
            faces[face_indices[4][0]]=children_left->faces[face_indices[5][1]];faces[face_indices[6][0]]=children_left->faces[face_indices[7][1]];}
        else{
            faces[face_indices[0][0]]=++number_of_faces;faces[face_indices[2][0]]=++number_of_faces;faces[face_indices[4][0]]=++number_of_faces;faces[face_indices[6][0]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],0));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[2],0));
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[4],0));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[6],0));}}
        if(children_right){
            faces[face_indices[1][1]]=children_right->faces[face_indices[0][0]];faces[face_indices[3][1]]=children_right->faces[face_indices[2][0]];
            faces[face_indices[5][1]]=children_right->faces[face_indices[4][0]];faces[face_indices[7][1]]=children_right->faces[face_indices[6][0]];}
        else{
            faces[face_indices[1][1]]=++number_of_faces;faces[face_indices[3][1]]=++number_of_faces;faces[face_indices[5][1]]=++number_of_faces;faces[face_indices[7][1]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[1],1));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[3],1));
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[5],1));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[7],1));}}
        if(children_bottom){
            faces[face_indices[0][2]]=children_bottom->faces[face_indices[2][3]];faces[face_indices[1][2]]=children_bottom->faces[face_indices[3][3]];
            faces[face_indices[4][2]]=children_bottom->faces[face_indices[6][3]];faces[face_indices[5][2]]=children_bottom->faces[face_indices[7][3]];}
        else{
            faces[face_indices[0][2]]=++number_of_faces;faces[face_indices[1][2]]=++number_of_faces;faces[face_indices[4][2]]=++number_of_faces;faces[face_indices[5][2]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],2));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[1],2));
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[4],2));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[5],2));}}
        if(children_top){
            faces[face_indices[2][3]]=children_top->faces[face_indices[0][2]];faces[face_indices[3][3]]=children_top->faces[face_indices[1][2]];
            faces[face_indices[6][3]]=children_top->faces[face_indices[4][2]];faces[face_indices[7][3]]=children_top->faces[face_indices[5][2]];}
        else{
            faces[face_indices[2][3]]=++number_of_faces;faces[face_indices[3][3]]=++number_of_faces;faces[face_indices[6][3]]=++number_of_faces;faces[face_indices[7][3]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[2],3));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[3],3));
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[6],3));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[7],3));}}
        if(children_front){
            faces[face_indices[0][4]]=children_front->faces[face_indices[4][5]];faces[face_indices[1][4]]=children_front->faces[face_indices[5][5]];
            faces[face_indices[2][4]]=children_front->faces[face_indices[6][5]];faces[face_indices[3][4]]=children_front->faces[face_indices[7][5]];}
        else{
            faces[face_indices[0][4]]=++number_of_faces;faces[face_indices[1][4]]=++number_of_faces;faces[face_indices[2][4]]=++number_of_faces;faces[face_indices[3][4]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[0],4));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[1],4));
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[4],4));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[5],4));}}
        if(children_back){
            faces[face_indices[4][5]]=children_back->faces[face_indices[0][4]];faces[face_indices[5][5]]=children_back->faces[face_indices[1][4]];
            faces[face_indices[6][5]]=children_back->faces[face_indices[2][4]];faces[face_indices[7][5]]=children_back->faces[face_indices[3][4]];}
        else{
            faces[face_indices[4][5]]=++number_of_faces;faces[face_indices[5][5]]=++number_of_faces;faces[face_indices[6][5]]=++number_of_faces;faces[face_indices[7][5]]=++number_of_faces;
            if(new_faces){
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[4],5));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[5],5));
                new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[6],5));new_faces->Append(PAIR<OCTREE_CELL<T>*,int>(&children[7],5));}}}
}
//#####################################################################
// Function Edge_Neighbor_Check
//#####################################################################
template<class T> void OCTREE_CHILDREN<T>::
Edge_Neighbor_Check(int& number_of_nodes,const int x_offset,const int y_offset,const int z_offset,const int node_index,const int neighbor_node_index,const int child_number,const int node_number,
                    ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes,const OCTREE_GRID<T>* grid)
{
    if(nodes[node_index] == 0){
        OCTREE_CELL<T>* neighbor=parent->Get_Neighbor(x_offset,y_offset,z_offset,grid);
        OCTREE_CHILDREN<T>* edge_neighbor=(neighbor==0)?0:neighbor->children;
        if(edge_neighbor) nodes[node_index]=edge_neighbor->nodes[neighbor_node_index];
        else{
            nodes[node_index]=++number_of_nodes;
            if(new_nodes) new_nodes->Append(PAIR<OCTREE_CELL<T>*,int>(&children[child_number],node_number));}}
}
//#####################################################################
// Function Storage_Requirement
//#####################################################################
template<class T> int OCTREE_CHILDREN<T>::
Storage_Requirement() const
{
    int total=sizeof(OCTREE_CELL<T>*)+2*sizeof(int*)*sizeof(VECTOR<T,3>)+sizeof(int);
    if(nodes)total+=27*sizeof(int);
    if(faces)total+=36*sizeof(int);
    for(int i=0;i<8;i++)total+=children[i].Storage_Requirement();
    return total;
}
//#####################################################################
template class OCTREE_CHILDREN<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OCTREE_CHILDREN<double>;
#endif
#endif
