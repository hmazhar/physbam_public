//#####################################################################
// Copyright 2004, Zhaosheng Bao, Ron Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRON_HIERARCHY<T>::
TETRAHEDRON_HIERARCHY(TETRAHEDRON_MESH& tetrahedron_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,const bool update_boxes)
    :tetrahedron_mesh(tetrahedron_mesh_input),particles(particles_input),tetrahedron_list(0)
{  
    if(tetrahedron_mesh.elements.m){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
} 
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRON_HIERARCHY<T>::
TETRAHEDRON_HIERARCHY(TETRAHEDRON_MESH& tetrahedron_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,ARRAY<TETRAHEDRON<T> >& tetrahedron_list_input,const bool update_boxes)
    :tetrahedron_mesh(tetrahedron_mesh_input),particles(particles_input),tetrahedron_list(&tetrahedron_list_input)
{  
    if(tetrahedron_mesh.elements.m){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TETRAHEDRON_HIERARCHY<T>::
~TETRAHEDRON_HIERARCHY()
{}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class T> void TETRAHEDRON_HIERARCHY<T>::
Initialize_Hierarchy_Using_KD_Tree()
{
    KD_TREE<TV> kd_tree(false);
    ARRAY<TV> centroids(tetrahedron_mesh.elements.m);
    for(int t=1;t<=tetrahedron_mesh.elements.m;t++){
        int i,j,k,l;tetrahedron_mesh.elements(t).Get(i,j,k,l);centroids(t)=TETRAHEDRON<T>::Center(particles.X(i),particles.X(j),particles.X(k),particles.X(l));}
    kd_tree.Create_Left_Balanced_KD_Tree(centroids);
    leaves=tetrahedron_mesh.elements.m;parents.Resize(leaves);children.Remove_All();
    root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);assert(root==2*leaves-1);
    box_hierarchy.Resize(root);box_radius.Resize(root);
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class T> void TETRAHEDRON_HIERARCHY<T>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> X) 
{
    for(int k=1;k<=leaves;k++){const VECTOR<int,4>& nodes=tetrahedron_mesh.elements(k);
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(X.Subset(nodes));}
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class T> void TETRAHEDRON_HIERARCHY<T>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> start_X,ARRAY_VIEW<const TV> end_X)
{
    for(int k=1;k<=leaves;k++){const VECTOR<int,4>& nodes=tetrahedron_mesh.elements(k);
        bounding_boxes(k)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(start_X.Subset(nodes)),RANGE<TV>::Bounding_Box(end_X.Subset(nodes)));}
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class T> void TETRAHEDRON_HIERARCHY<T>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const FRAME<TV>& start_frame,const FRAME<TV>& end_frame)
{
    for(int k=1;k<=leaves;k++){const VECTOR<int,4>& nodes=tetrahedron_mesh.elements(k);
        bounding_boxes(k)=RANGE<TV>::Combine(
            RANGE<TV>::Bounding_Box(start_frame*particles.X(nodes[1]),start_frame*particles.X(nodes[2]),start_frame*particles.X(nodes[3]),start_frame*particles.X(nodes[4])),
            RANGE<TV>::Bounding_Box(end_frame*particles.X(nodes[1]),end_frame*particles.X(nodes[2]),end_frame*particles.X(nodes[3]),end_frame*particles.X(nodes[4])));}
}
//#####################################################################
// Function Calculate_Bounding_Box_Radii
//#####################################################################
// at the boxes center, but may be a tighter bound on the tetrahedron than the box
template<class T> void TETRAHEDRON_HIERARCHY<T>::
Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) 
{
    for(int k=1;k<=leaves;k++){
        TV center=bounding_boxes(k).Center();int node1,node2,node3,node4;tetrahedron_mesh.elements(k).Get(node1,node2,node3,node4);
        radius(k)=sqrt(max((particles.X(node1)-center).Magnitude_Squared(),(particles.X(node2)-center).Magnitude_Squared(),(particles.X(node3)-center).Magnitude_Squared(),
            (particles.X(node4)-center).Magnitude_Squared()));}
}
//#####################################################################
template class TETRAHEDRON_HIERARCHY<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TETRAHEDRON_HIERARCHY<double>;
#endif
