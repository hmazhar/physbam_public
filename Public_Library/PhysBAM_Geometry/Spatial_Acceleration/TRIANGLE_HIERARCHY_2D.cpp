//#####################################################################
// Copyright 2002-2004, Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class T> void TRIANGLE_HIERARCHY_2D<T>::
Initialize_Hierarchy_Using_KD_Tree()
{
    KD_TREE<TV> kd_tree(false);
    ARRAY<TV> centroids(triangle_mesh.elements.m);
    for(int t=1;t<=triangle_mesh.elements.m;t++){int i,j,k;triangle_mesh.elements(t).Get(i,j,k);centroids(t)=TRIANGLE_2D<T>::Center(particles.X(i),particles.X(j),particles.X(k));}
    kd_tree.Create_Left_Balanced_KD_Tree(centroids);
    leaves=triangle_mesh.elements.m;parents.Resize(leaves);children.Remove_All();
    root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);assert(root==2*leaves-1);
    box_hierarchy.Resize(root);box_radius.Resize(root);
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class T> void TRIANGLE_HIERARCHY_2D<T>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> X) 
{
    for(int k=1;k<=leaves;k++){const VECTOR<int,3>& nodes=triangle_mesh.elements(k);
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(X.Subset(nodes));}
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class T> void TRIANGLE_HIERARCHY_2D<T>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> start_X,ARRAY_VIEW<const TV> end_X)
{
    for(int k=1;k<=leaves;k++){const VECTOR<int,3>& nodes=triangle_mesh.elements(k);
        bounding_boxes(k)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(start_X.Subset(nodes)),RANGE<TV>::Bounding_Box(end_X.Subset(nodes)));}
}
//#####################################################################
// Function Calculate_Bounding_Box_Radii
//#####################################################################
// at the boxes center, but may be a tighter bound on the triangle than the box
template<class T> void TRIANGLE_HIERARCHY_2D<T>::
Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) 
{
    for(int k=1;k<=leaves;k++){
        TV center=bounding_boxes(k).Center();int node1,node2,node3;triangle_mesh.elements(k).Get(node1,node2,node3);
        radius(k)=sqrt(max((particles.X(node1)-center).Magnitude_Squared(),(particles.X(node2)-center).Magnitude_Squared(),(particles.X(node3)-center).Magnitude_Squared()));}
}
//#####################################################################
template class TRIANGLE_HIERARCHY_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_HIERARCHY_2D<double>;
#endif
