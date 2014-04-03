//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Geometry/Red_Green/RED_TETRAHEDRON.h>
using namespace PhysBAM;
//#####################################################################
template<class T> const TETRAHEDRAL_GROUP<T> RED_CHILDREN_3D<T>::relative_orientations[8]={0,0,0,0,TETRAHEDRAL_GROUP<T>::a,TETRAHEDRAL_GROUP<T>::a2,TETRAHEDRAL_GROUP<T>::b,TETRAHEDRAL_GROUP<T>::d2};
template<class T> const int RED_CHILDREN_3D<T>::node_indices[8][4]={{0,4,6,7},{4,1,9,8},{6,9,2,5},{7,8,5,3},{6,7,4,5},{4,5,7,8},{8,9,4,5},{4,5,9,6}};
//#####################################################################
// Function Create_Node_Compaction_Mapping_Helper
//#####################################################################
template<class T> void RED_CHILDREN_3D<T>::
Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count)
{
    for(int i=0;i<10;i++)if(nodes[i]){
        if(!mapping_array(nodes[i]))mapping_array(nodes[i])=++node_count;
        nodes[i]=mapping_array(nodes[i]);}
    for(int i=0;i<8;i++){
        if(children[i].Has_Red_Children()) children[i].red_children->Create_Node_Compaction_Mapping_Helper(mapping_array,node_count);
        else if(children[i].Has_Green_Children()) children[i].green_children->Create_Node_Compaction_Mapping_Helper(mapping_array,node_count);}
}
//#####################################################################
// Function Initialize_Pseudo_Root_Tetrahedron
//#####################################################################
template<class T> void RED_CHILDREN_3D<T>::
Initialize_Pseudo_Root_Tetrahedron(int& number_of_cells,int node1,int node2,int node3,int node4,const VECTOR<T,3>& center_input,const T childrens_size_input,const TETRAHEDRAL_GROUP<T>& orientation_input)
{
    parent=0;childrens_depth=1;center=center_input;childrens_size=childrens_size_input;orientation=orientation_input;
    children[0].Initialize(this,number_of_cells,0);center+=center-children[0].Center();
    for(int i=1;i<8;i++)children[i].Initialize(this);
    for(int i=0;i<10;i++)nodes[i]=0;
    nodes[node_indices[0][0]]=node1;nodes[node_indices[0][1]]=node2;nodes[node_indices[0][2]]=node3;nodes[node_indices[0][3]]=node4;
}
//#####################################################################
// Function Initialize_Non_Root_Tetrahedron
//#####################################################################
template<class T> void RED_CHILDREN_3D<T>::
Initialize_Non_Root_Tetrahedron(RED_TETRAHEDRON<T>* parent_input,int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid)
{
    parent=parent_input;childrens_depth=parent->owner->childrens_depth+1;
    center=parent->Center();childrens_size=(T).5*parent->owner->childrens_size;orientation=parent->Orientation();
    for(int i=0;i<8;i++)children[i].Initialize(this,number_of_cells,i);
    for(int i=0;i<4;i++)Corner(i)=parent->Node(i); // corner nodes from the parent tetrahedron
    for(int i=0;i<6;i++){ // add new midpoint nodes
        int neighbor_midpoint=parent->Get_Neighbor_Midpoint(i,grid);
        Midpoint(i)=neighbor_midpoint?neighbor_midpoint:++number_of_nodes;}
}
//#####################################################################
template class RED_CHILDREN_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RED_CHILDREN_3D<double>;
#endif
#endif
