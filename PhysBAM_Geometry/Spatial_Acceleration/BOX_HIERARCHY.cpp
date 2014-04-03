//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Joseph Teran, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_LINE_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_POINT_SIMPLEX_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOX_HIERARCHY<TV>::
BOX_HIERARCHY()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BOX_HIERARCHY<TV>::
~BOX_HIERARCHY()
{}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Clean_Memory()
{
    leaves=root=0;
    parents.Clean_Memory();children.Clean_Memory();box_hierarchy.Clean_Memory();traversal_stack.Clean_Memory();dual_traversal_stack.Clean_Memory();box_radius.Clean_Memory();
}
//#####################################################################
// Function Set_Leaf_Boxes
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Set_Leaf_Boxes(const ARRAY<RANGE<TV> >& boxes,const bool reinitialize)
{
    if(reinitialize){box_hierarchy=boxes;Initialize_Hierarchy_Using_KD_Tree();}
    else box_hierarchy.Prefix(leaves)=boxes;
    Update_Nonleaf_Boxes();
}
//#####################################################################
// Function Thicken_Leaf_Boxes
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Thicken_Leaf_Boxes(const T extra_thickness)
{
    for(int k=1;k<=leaves;k++) box_hierarchy(k).Change_Size(extra_thickness);
}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Initialize_Hierarchy_Using_KD_Tree()
{
    KD_TREE<TV> kd_tree(false);
    ARRAY<TV> centroids(box_hierarchy.m);
    for(int l=1;l<=box_hierarchy.m;l++)centroids(l)=box_hierarchy(l).Center();
    kd_tree.Create_Left_Balanced_KD_Tree(centroids);
    leaves=box_hierarchy.m;parents.Resize(leaves);children.Remove_All();
    root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);assert(root==2*leaves-1);
    box_hierarchy.Resize(root);
}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree_Helper
//#####################################################################
template<class TV> int BOX_HIERARCHY<TV>::
Initialize_Hierarchy_Using_KD_Tree_Helper(KD_TREE_NODE<T>* node)
{
    if(!node->left&&!node->right)return node->node_index;
    int left_child=Initialize_Hierarchy_Using_KD_Tree_Helper(node->left);
    if(!node->right)return left_child;
    int right_child=Initialize_Hierarchy_Using_KD_Tree_Helper(node->right);
    children.Append(VECTOR<int,2>(left_child,right_child));parents.Append(0);
    return parents(left_child)=parents(right_child)=children.m+leaves;
}
//#####################################################################
// Function Update_Nonleaf_Boxes
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Update_Nonleaf_Boxes()
{
    for(int k=leaves+1;k<=box_hierarchy.m;k++)
        box_hierarchy(k)=RANGE<TV>::Combine(box_hierarchy(children(k-leaves)(1)),box_hierarchy(children(k-leaves)(2)));
}
//#####################################################################
// Function Update_Modified_Nonleaf_Boxes
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Update_Modified_Nonleaf_Boxes(ARRAY<bool>& modified)
{
    for(int k=leaves+1;k<=box_hierarchy.m;k++){
        const VECTOR<int,2>& child=children(k-leaves);
        modified(k)=modified(child[1]) || modified(child[2]);
        if(modified(k)) box_hierarchy(k)=RANGE<TV>::Combine(box_hierarchy(child[1]),box_hierarchy(child[2]));}
}
//#####################################################################
// Function Calculate_Bounding_Box_Radii
//#####################################################################
// at the boxes center
template<class TV> void BOX_HIERARCHY<TV>::
Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) 
{
    for(int k=1;k<=leaves;k++) radius(k)=(T).5*bounding_boxes(k).Edge_Lengths().Magnitude();
}
//#####################################################################
// Function Update_Nonleaf_Box_Radii
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Update_Nonleaf_Box_Radii()
{
    for(int k=leaves+1;k<=box_radius.m;k++){
        int box1,box2;children(k-leaves).Get(box1,box2);TV center=box_hierarchy(k).Center();
        box_radius(k)=max((box_hierarchy(box1).Center()-center).Magnitude()+box_radius(box1),(box_hierarchy(box2).Center()-center).Magnitude()+box_radius(box2));}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_THICKNESS> void BOX_HIERARCHY<TV>::
Intersection_List(const int box,const TV& point,ARRAY<int>& intersection_list,const T_THICKNESS thickness_over_two) const
{
    if(!IS_SAME<T_THICKNESS,ZERO>::value && !thickness_over_two){Intersection_List(box,point,intersection_list,ZERO());return;}
    if(box==0) return;
    traversal_stack.Remove_All();traversal_stack.Push(box);
    while(!traversal_stack.Empty()){int current=traversal_stack.Pop();
        if(box_hierarchy(current).Outside(point,thickness_over_two)) continue;
        if(Leaf(current)) intersection_list.Append(current);else{int box1,box2;children(current-leaves).Get(box1,box2);traversal_stack.Push(box1);traversal_stack.Push(box2);}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_THICKNESS> void BOX_HIERARCHY<TV>::
Intersection_List(const int box,const RANGE<TV>& test_box,ARRAY<int>& intersection_list,const T_THICKNESS thickness_over_two) const
{
    if(!IS_SAME<T_THICKNESS,ZERO>::value && !thickness_over_two){Intersection_List(box,test_box,intersection_list,ZERO());return;}
    if(box==0) return;
    traversal_stack.Remove_All();traversal_stack.Push(box);
    while(!traversal_stack.Empty()){int current=traversal_stack.Pop();
        if(!test_box.Intersection(box_hierarchy(current),thickness_over_two)) continue;
        if(Leaf(current)) intersection_list.Append(current);else{int box1,box2;children(current-leaves).Get(box1,box2);traversal_stack.Push(box1);traversal_stack.Push(box2);}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Intersection_List(const int box,const ORIENTED_BOX<TV>& test_box,ARRAY<int>& intersection_list) const
{
    traversal_stack.Remove_All();traversal_stack.Push(box);
    while(!traversal_stack.Empty()){int current=traversal_stack.Pop();
        if(!test_box.Intersection(box_hierarchy(current))) continue;
        if(Leaf(current)) intersection_list.Append(current);else{int box1,box2;children(current-leaves).Get(box1,box2);traversal_stack.Push(box1);traversal_stack.Push(box2);}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Intersection_List(const int box,const T_HYPERPLANE& test_plane,ARRAY<int>& intersection_list,const T thickness_over_two) const
{
    traversal_stack.Remove_All();traversal_stack.Push(box);
    while(!traversal_stack.Empty()){int current=traversal_stack.Pop();
        if(!INTERSECTION::Intersects(box_hierarchy(current),test_plane,thickness_over_two)) continue;
        if(Leaf(current)) intersection_list.Append(current);else{int box1,box2;children(current-leaves).Get(box1,box2);traversal_stack.Push(box1);traversal_stack.Push(box2);}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> void BOX_HIERARCHY<TV>::
Intersection_List(const int box,const IMPLICIT_OBJECT<TV>& implicit_object,const MATRIX<T,TV::dimension>& rotation,const TV& translation,ARRAY<int>& intersection_list,const T contour_value) const
{
    PHYSBAM_ASSERT(box_radius.m==box_hierarchy.m);
    traversal_stack.Remove_All();traversal_stack.Push(box);
    while(!traversal_stack.Empty()){int current=traversal_stack.Pop();
        TV location=rotation*box_hierarchy(current).Center()+translation; // location in space of the implicit object
        if(implicit_object.Lazy_Outside_Extended_Levelset(location,box_radius(current)+contour_value)) continue;
        if(Leaf(current)) intersection_list.Append(current);else{int box1,box2;children(current-leaves).Get(box1,box2);traversal_stack.Push(box1);traversal_stack.Push(box2);}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class BOX_HIERARCHY<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif

