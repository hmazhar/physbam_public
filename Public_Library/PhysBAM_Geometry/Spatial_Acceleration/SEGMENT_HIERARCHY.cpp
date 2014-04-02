//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SEGMENT_HIERARCHY<TV>::
SEGMENT_HIERARCHY(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,const bool update_boxes)
    :segment_mesh(segment_mesh_input),particles(particles_input),segment_list(0)
{
    if(segment_mesh.elements.m){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SEGMENT_HIERARCHY<TV>::
SEGMENT_HIERARCHY(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,ARRAY<T_SEGMENT>& segment_list_input,const bool update_boxes)
    :segment_mesh(segment_mesh_input),particles(particles_input),segment_list(&segment_list_input)
{
    if(segment_mesh.elements.m){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SEGMENT_HIERARCHY<TV>::
~SEGMENT_HIERARCHY()
{}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class TV> void SEGMENT_HIERARCHY<TV>::
Initialize_Hierarchy_Using_KD_Tree()
{
    KD_TREE<TV> kd_tree(false);
    ARRAY<TV> centroids(segment_mesh.elements.m);
    for(int t=1;t<=segment_mesh.elements.m;t++){
        int i,j;segment_mesh.elements(t).Get(i,j);centroids(t)=(T).5*(particles.X(i)+particles.X(j));}
    kd_tree.Create_Left_Balanced_KD_Tree(centroids);
    leaves=segment_mesh.elements.m;parents.Resize(leaves);children.Remove_All();
    root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);assert(root==2*leaves-1);
    box_hierarchy.Resize(root);box_radius.Resize(root);
}
//#####################################################################
// Function Calculate_Bounding_Boxes_Helper
//#####################################################################
template<class TV> template<class T_ARRAY_TV> void SEGMENT_HIERARCHY<TV>::
Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV X)
{
    STATIC_ASSERT((IS_SAME<TV,typename T_ARRAY_TV::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY_TV>::value));
    for(int k=1;k<=leaves;k++){
        int node1,node2;segment_mesh.elements(k).Get(node1,node2);
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(X(node1),X(node2));}
}
//#####################################################################
// Function Calculate_Bounding_Boxes_Helper
//#####################################################################
template<class TV> template<class T_ARRAY_TV> void SEGMENT_HIERARCHY<TV>::
Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV start_X,T_ARRAY_TV end_X)
{
    STATIC_ASSERT((IS_SAME<TV,typename T_ARRAY_TV::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY_TV>::value));
    for(int k=1;k<=leaves;k++){
        int node1,node2;segment_mesh.elements(k).Get(node1,node2);
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(start_X(node1),start_X(node2),end_X(node1),end_X(node2));}
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class TV> void SEGMENT_HIERARCHY<TV>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const FRAME<TV>& start_frame,const FRAME<TV>& end_frame)
{
    for(int k=1;k<=leaves;k++){
        int node1,node2;segment_mesh.elements(k).Get(node1,node2);
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(start_frame*particles.X(node1),start_frame*particles.X(node2),end_frame*particles.X(node1),end_frame*particles.X(node2));}
}
//#####################################################################
// Function Calculate_Bounding_Box_Radii
//#####################################################################
// at the boxes center, but may be a tighter bound on the segment than the box
template<class TV> void SEGMENT_HIERARCHY<TV>::
Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius)
{
    for(int k=1;k<=leaves;k++){
        TV center=bounding_boxes(k).Center();int node1,node2;segment_mesh.elements(k).Get(node1,node2);
        radius(k)=sqrt(max((particles.X(node1)-center).Magnitude_Squared(),(particles.X(node2)-center).Magnitude_Squared()));}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class SEGMENT_HIERARCHY<VECTOR<T,d> >; \
    template void SEGMENT_HIERARCHY<VECTOR<T,d> >::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<VECTOR<T,d> > >&,ARRAY_VIEW<const VECTOR<T,d> >); \
    template void SEGMENT_HIERARCHY<VECTOR<T,d> >::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<VECTOR<T,d> > >&,INDIRECT_ARRAY<ARRAY_VIEW<const VECTOR<T,d> > >);
INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
template void SEGMENT_HIERARCHY<VECTOR<float,2> >::Calculate_Bounding_Boxes_Helper<ARRAY_VIEW<VECTOR<float,2> const,int> >(ARRAY<RANGE<VECTOR<float,2> >,int>&,
    ARRAY_VIEW<VECTOR<float,2> const,int>,ARRAY_VIEW<VECTOR<float,2> const,int>);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
template void SEGMENT_HIERARCHY<VECTOR<double,2> >::Calculate_Bounding_Boxes_Helper<ARRAY_VIEW<VECTOR<double,2> const,int> >(ARRAY<RANGE<VECTOR<double,2> >,int>&,
    ARRAY_VIEW<VECTOR<double,2> const,int>,ARRAY_VIEW<VECTOR<double,2> const,int>);
#endif
