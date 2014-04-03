//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T_ARRAY> PARTICLE_HIERARCHY<TV,T_ARRAY>::
PARTICLE_HIERARCHY(const T_ARRAY& X_input,const bool update_boxes,const int particles_per_group_input)
    :X(X_input),particles_per_group(particles_per_group_input)
{
    if(X.Size()){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T_ARRAY> PARTICLE_HIERARCHY<TV,T_ARRAY>::
~PARTICLE_HIERARCHY()
{}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Initialize_Hierarchy_Using_KD_Tree()
{
    particles_in_group.Clean_Memory();
    KD_TREE<TV> kd_tree(false);
    ARRAY<TV> X_copy(X);
    if(particles_per_group){
        kd_tree.Create_Left_Balanced_KD_Tree_With_Grouping(X_copy,particles_in_group,particles_per_group);
        leaves=particles_in_group.m;}
    else{
        kd_tree.Create_Left_Balanced_KD_Tree(X_copy);
        leaves=X.Size();}
    parents.Resize(leaves);
    children.Remove_All();root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);
    assert(root==2*leaves-1);box_hierarchy.Resize(root);box_radius.Resize(root);
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes) const
{
    if(particles_per_group) for(int k=1;k<=leaves;++k){
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(X.Subset(particles_in_group(k)));}
    else for(int k=1;k<=leaves;++k) bounding_boxes(k)=RANGE<TV>(X(k),X(k));
}
//#####################################################################
// Function Calculate_Bounding_Box_Radii
//#####################################################################
// at the boxes center, but may be a tighter bound on the triangle than the box
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius)
{
    if(particles_per_group) for(int k=1;k<=leaves;k++){
        TV center=bounding_boxes(k).Center();T max_radius_squared=0;
        for(int i=1;i<=particles_in_group(k).m;i++) max_radius_squared=max(max_radius_squared,(X(particles_in_group(k)(i))-center).Magnitude_Squared());
        radius(k)=sqrt(max_radius_squared);}
    else for(int k=1;k<=leaves;k++){TV center=bounding_boxes(k).Center();radius(k)=sqrt((X(k)-center).Magnitude_Squared());}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class PARTICLE_HIERARCHY<VECTOR<T,d> >; \
    template class PARTICLE_HIERARCHY<VECTOR<T,d>,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,d>,int> > >;
INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif
