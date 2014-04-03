//#####################################################################
// Copyright 2004-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KD_TREE
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE_NODE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <algorithm>
using namespace PhysBAM;
//#####################################################################
// Function Inside_Node
//#####################################################################
template<class TV> KD_TREE_NODE<typename TV::SCALAR>* KD_TREE<TV>::
Leaf_Node(const TV& point) const
{
    return Leaf_Node_Helper(point,root_node);
}
//#####################################################################
// Function Inside_Node_Helper
//#####################################################################
template<class TV> KD_TREE_NODE<typename TV::SCALAR>* KD_TREE<TV>::
Leaf_Node_Helper(const TV& point,KD_TREE_NODE<T>* node) const
{
    if(!node->split_axis) return node;
    else{
        if(point[node->split_axis]<node->split_value) return Leaf_Node_Helper(point,node->left);
        else if(point[node->split_axis]>node->split_value) return Leaf_Node_Helper(point,node->right);
        else return node;}
}
//#####################################################################
// Function Create_KD_Tree
//#####################################################################
template<class TV> void KD_TREE<TV>::
Create_KD_Tree(ARRAY_VIEW<const TV> points)
{
    PHYSBAM_ASSERT(points.Size(),"trying to create empty kd-tree");
    if(store_values_on_internal_nodes){
        pool.Delete_All();root_node=pool.New();
        ARRAY<int> permutation_array(IDENTITY_ARRAY<>(points.Size()));
        RANGE<TV> box=RANGE<TV>::Bounding_Box(points);
        Balance_Sub_KD_Tree_Using_Internal_Nodes(root_node,1,permutation_array.m,points,permutation_array,box);}
    else{
        pool.Delete_All();root_node=pool.New();
        ARRAY<int> permutation_array(IDENTITY_ARRAY<>(points.Size()));
        RANGE<TV> box=RANGE<TV>::Bounding_Box(points);
        Sub_KD_Tree_Using_Leaf_Nodes(root_node,1,permutation_array.m,points,permutation_array,box);}
}
//#####################################################################
// Function Create_Left_Balanced_KD_Tree
//#####################################################################
template<class TV> void KD_TREE<TV>::
Create_Left_Balanced_KD_Tree(ARRAY_VIEW<const TV> points_to_balance)
{
    PHYSBAM_ASSERT(points_to_balance.Size(),"trying to create empty kd-tree");
    pool.Delete_All();root_node=pool.New();
    ARRAY<int> permutation_array(IDENTITY_ARRAY<>(points_to_balance.Size()));
    RANGE<TV> box=RANGE<TV>::Bounding_Box(points_to_balance);
    if(store_values_on_internal_nodes) Balance_Sub_KD_Tree_Using_Internal_Nodes(root_node,1,permutation_array.m,points_to_balance,permutation_array,box);
    else Balance_Sub_KD_Tree_Using_Leaf_Nodes(root_node,1,permutation_array.m,points_to_balance,permutation_array,box);
}
//#####################################################################
// Function Create_Left_Balanced_KD_Tree_With_Grouping
//#####################################################################
// Leaf index refers to a group number rather than an individual point.  points_in_group lists group elements.
template<class TV> void KD_TREE<TV>::
Create_Left_Balanced_KD_Tree_With_Grouping(ARRAY_VIEW<const TV> points_to_balance,ARRAY<ARRAY<int> >& points_in_group,const int max_points_in_group)
{
    PHYSBAM_ASSERT(points_to_balance.Size(),"trying to create empty kd-tree");
    pool.Delete_All();root_node=pool.New();
    ARRAY<int> permutation_array((IDENTITY_ARRAY<>(points_to_balance.Size())));
    ARRAY<int> point_group(points_to_balance.Size());int current_group_index=0;
    RANGE<TV> box=RANGE<TV>::Bounding_Box(points_to_balance);
    PHYSBAM_ASSERT(!store_values_on_internal_nodes);
    Balance_Sub_KD_Tree_Using_Leaf_Nodes_With_Grouping(root_node,1,points_to_balance.Size(),points_to_balance,permutation_array,box,point_group,current_group_index,max_points_in_group);
    int number_of_groups=ARRAYS_COMPUTATIONS::Max(point_group);
    ARRAY<int> group_size(number_of_groups);for(int i=1;i<=point_group.m;i++) group_size(point_group(i))++;
    points_in_group.Resize(number_of_groups);for(int i=1;i<=points_in_group.m;i++) points_in_group(i).Resize(group_size(i));
    for(int i=1;i<=points_to_balance.Size();i++) points_in_group(point_group(i))(group_size(point_group(i))--)=i;
}
//#####################################################################
// Function Balance_Sub_KD_Tree_Using_Internal_Nodes
//#####################################################################
template<class TV> void KD_TREE<TV>::
Balance_Sub_KD_Tree_Using_Internal_Nodes(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,RANGE<TV>& box)
{
    if(last_index==first_index){cell->split_axis=0;cell->node_index=permutation_array(first_index);return;}
    int partition_index=Choose_Partition_Index_Using_Internal_Nodes(first_index,last_index);
    int axis=Choose_Partition_Axis(box.Edge_Lengths());
    Median_Split(partition_index,first_index,last_index,points,permutation_array,axis);
    cell->split_axis=axis;cell->split_value=points(permutation_array(partition_index))[axis];cell->node_index=permutation_array(partition_index);
    if(partition_index>first_index){
        cell->left=pool.New();
        RANGE<TV> left_box(box);left_box.max_corner[axis]=cell->split_value;
        Balance_Sub_KD_Tree_Using_Internal_Nodes(cell->left,first_index,partition_index-1,points,permutation_array,left_box);}
    if(partition_index<last_index){
        cell->right=pool.New();
        RANGE<TV> right_box(box);right_box.min_corner[axis]=cell->split_value;
        Balance_Sub_KD_Tree_Using_Internal_Nodes(cell->right,partition_index+1,last_index,points,permutation_array,right_box);}
}
//#####################################################################
// Function Balance_Sub_KD_Tree_Using_Leaf_Nodes
//#####################################################################
template<class TV> void KD_TREE<TV>::
Balance_Sub_KD_Tree_Using_Leaf_Nodes(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,RANGE<TV>& box)
{
    if(last_index==first_index){cell->split_axis=0;cell->node_index=permutation_array(first_index);return;}
    int partition_index=Choose_Partition_Index_Using_Leaf_Nodes(first_index,last_index);
    cell->split_axis=Choose_Partition_Axis(box.Edge_Lengths());
    Median_Split(partition_index,first_index,last_index,points,permutation_array,cell->split_axis);
    cell->split_value=points(permutation_array(partition_index))[cell->split_axis];
    assert(partition_index>first_index);
    cell->left=pool.New();
    RANGE<TV> left_box(box);left_box.max_corner[cell->split_axis]=cell->split_value;
    Balance_Sub_KD_Tree_Using_Leaf_Nodes(cell->left,first_index,partition_index-1,points,permutation_array,left_box);
    if(partition_index<=last_index){
        cell->right=pool.New();
        RANGE<TV> right_box(box);right_box.min_corner[cell->split_axis]=cell->split_value;
        Balance_Sub_KD_Tree_Using_Leaf_Nodes(cell->right,partition_index,last_index,points,permutation_array,right_box);}
}
//#####################################################################
// Function Balance_Sub_KD_Tree_Using_Leaf_Nodes_With_Grouping
//#####################################################################
template<class TV> void KD_TREE<TV>::
Balance_Sub_KD_Tree_Using_Leaf_Nodes_With_Grouping(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,
    RANGE<TV>& box,ARRAY<int>& point_group,int& current_group_index,const int max_points_in_group)
{
    if(last_index-first_index+1<=max_points_in_group){
        cell->split_axis=0;cell->node_index=++current_group_index;
        for(int i=first_index;i<=last_index;i++) point_group(permutation_array(i))=current_group_index;return;}
    int partition_index=Choose_Partition_Index_Using_Leaf_Nodes(first_index,last_index);
    cell->split_axis=Choose_Partition_Axis(box.Edge_Lengths());
    Median_Split(partition_index,first_index,last_index,points,permutation_array,cell->split_axis);
    cell->split_value=points(permutation_array(partition_index))[cell->split_axis];
    assert(partition_index>first_index);
    cell->left=pool.New();
    RANGE<TV> left_box(box);left_box.max_corner[cell->split_axis]=cell->split_value;
    Balance_Sub_KD_Tree_Using_Leaf_Nodes_With_Grouping(cell->left,first_index,partition_index-1,points,permutation_array,left_box,point_group,current_group_index,max_points_in_group);
    if(partition_index<=last_index){
        cell->right=pool.New();
        RANGE<TV> right_box(box);right_box.min_corner[cell->split_axis]=cell->split_value;
        Balance_Sub_KD_Tree_Using_Leaf_Nodes_With_Grouping(cell->right,partition_index,last_index,points,permutation_array,right_box,point_group,current_group_index,max_points_in_group);}
}
//#####################################################################
// Function Sub_KD_Tree_Using_Leaf_Nodes
//#####################################################################
template<class TV> void KD_TREE<TV>::
Sub_KD_Tree_Using_Leaf_Nodes(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,RANGE<TV>& box)
{
    if(last_index==first_index){cell->split_axis=0;cell->node_index=permutation_array(first_index);return;}
    int partition_index=(last_index+first_index+1)/2; // middle index right biased
    cell->split_axis=Choose_Partition_Axis(box.Edge_Lengths());
    Median_Split(partition_index,first_index,last_index,points,permutation_array,cell->split_axis);
    cell->split_value=points(permutation_array(partition_index))[cell->split_axis];
    assert(partition_index>first_index);
    cell->left=pool.New();
    RANGE<TV> left_box(box);left_box.max_corner[cell->split_axis]=cell->split_value;
    Sub_KD_Tree_Using_Leaf_Nodes(cell->left,first_index,partition_index-1,points,permutation_array,left_box);
    if(partition_index<=last_index){
        cell->right=pool.New();
        RANGE<TV> right_box(box);right_box.min_corner[cell->split_axis]=cell->split_value;
        Sub_KD_Tree_Using_Leaf_Nodes(cell->right,partition_index,last_index,points,permutation_array,right_box);}
}
//#####################################################################
// Function Median_Split
//#####################################################################
template<class TV> void KD_TREE<TV>::
Median_Split(const int partition_index,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,const int axis)
{
    ARRAY_VIEW<int> permutation_subset(permutation_array.Subset(first_index-1+IDENTITY_ARRAY<>(last_index-first_index+1)));
    ARRAY<T> values(points.Subset(permutation_subset).Project(axis)); // copy so that nth_element doesn't mess up original values
    std::nth_element(values.begin(),&values(partition_index-first_index+1),values.end());
    T split_value=values(partition_index-first_index+1);
    Partition_Helper_Less_Equal partition_helper_less_equal(&points,axis,split_value);
    Partition_Helper_Less partition_helper_less(&points,axis,split_value);
    int* middle=std::stable_partition(permutation_subset.begin(),permutation_subset.end(),partition_helper_less);
    std::stable_partition(middle,permutation_subset.end(),partition_helper_less_equal);
    assert(points(permutation_array(partition_index-1))[axis]<=split_value && split_value<=points(permutation_array(partition_index))[axis]);
}
//#####################################################################
// Function Find_Points_Within_Radius
//#####################################################################
template<class TV> void KD_TREE<TV>::
Find_Points_Within_Radius(const TV& location,const T max_distance_squared,ARRAY<int>& points_found,ARRAY<T>& distance_squared_of_points_found,
    ARRAY_VIEW<const TV> original_points) const
{
    if(!root_node) return;
    Find_Points_Within_Radius_Helper(root_node,location,max_distance_squared,points_found,distance_squared_of_points_found,original_points);
}
//#####################################################################
// Function Locate_Nearest_Neighbors
//#####################################################################
template<class TV> void KD_TREE<TV>::
Locate_Nearest_Neighbors(const TV& location,const T max_distance_squared,ARRAY<int>& points_found,ARRAY<T>& distance_squared_of_points_found,
    int& number_points_found,T& max_distance_squared_of_found_points,ARRAY_VIEW<const TV> original_points) const
{
    assert(points_found.m==distance_squared_of_points_found.m);assert(points_found.m>0);
    number_points_found=0;
    if(!root_node) return;
    T temp_max_distance=max_distance_squared;
    Locate_Nearest_Neighbors_Helper(root_node,location,temp_max_distance,number_points_found,points_found,distance_squared_of_points_found,original_points);
    if(number_points_found<=points_found.m) // otherwise max_distance_squared_of_found_points is computed in Locate_Nearest_Neighbors_Helper
        max_distance_squared_of_found_points=ARRAYS_COMPUTATIONS::Max(distance_squared_of_points_found.Prefix(number_points_found));
    else max_distance_squared_of_found_points=temp_max_distance;
    number_points_found=min(points_found.m,number_points_found);
}
//#####################################################################
// Function Find_Points_Within_Radius_Helper
//#####################################################################
template<class TV> void KD_TREE<TV>::
Find_Points_Within_Radius_Helper(const KD_TREE_NODE<T>* cell,const TV& location,const T max_distance_squared,ARRAY<int>& points_found,
    ARRAY<T>& distance_squared_of_points_found,ARRAY_VIEW<const TV> original_points) const
{
    if(cell->left || cell->right){
        T axis_distance=location[cell->split_axis]-cell->split_value;
        if(axis_distance>0){ // point belongs on right subtree
            if(cell->right) Find_Points_Within_Radius_Helper(cell->right,location,max_distance_squared,points_found,distance_squared_of_points_found,original_points);
            if(sqr(axis_distance)<max_distance_squared)
                if(cell->left) Find_Points_Within_Radius_Helper(cell->left,location,max_distance_squared,points_found,distance_squared_of_points_found,original_points);}
        else{ // point belongs on left subtree
            if(cell->left) Find_Points_Within_Radius_Helper(cell->left,location,max_distance_squared,points_found,distance_squared_of_points_found,original_points);
            if(sqr(axis_distance)<max_distance_squared)
                if(cell->right) Find_Points_Within_Radius_Helper(cell->right,location,max_distance_squared,points_found,distance_squared_of_points_found,original_points);}}

    T distance_squared_to_query_location=(original_points(cell->node_index)-location).Magnitude_Squared();
    if(distance_squared_to_query_location<=max_distance_squared){
        points_found.Append(cell->node_index);
        distance_squared_of_points_found.Append(distance_squared_to_query_location);}
}
//#####################################################################
// Function Locate_Nearest_Neighbors_Helper
//#####################################################################
template<class TV> void KD_TREE<TV>::
Locate_Nearest_Neighbors_Helper(const KD_TREE_NODE<T>* cell,const TV& location,T& max_distance_squared,int& number_of_points_found,ARRAY<int>& points_found,
    ARRAY<T>& distance_squared_of_points_found,ARRAY_VIEW<const TV> original_points) const
{
    if(cell->left || cell->right){
        T axis_distance=location[cell->split_axis]-cell->split_value;
        if(axis_distance>0){ // point belongs on right subtree
            if(cell->right) Locate_Nearest_Neighbors_Helper(cell->right,location,max_distance_squared,number_of_points_found,points_found,distance_squared_of_points_found,original_points);
            if(sqr(axis_distance)<max_distance_squared)
                if(cell->left) Locate_Nearest_Neighbors_Helper(cell->left,location,max_distance_squared,number_of_points_found,points_found,distance_squared_of_points_found,original_points);}
        else{ // point belongs on left subtree
            if(cell->left) Locate_Nearest_Neighbors_Helper(cell->left,location,max_distance_squared,number_of_points_found,points_found,distance_squared_of_points_found,original_points);
            if(sqr(axis_distance)<max_distance_squared)
                if(cell->right) Locate_Nearest_Neighbors_Helper(cell->right,location,max_distance_squared,number_of_points_found,points_found,distance_squared_of_points_found,original_points);}}

    T distance_squared_to_query_location=(original_points(cell->node_index)-location).Magnitude_Squared();
    if(distance_squared_to_query_location<=max_distance_squared){
        if(number_of_points_found<points_found.m){
            number_of_points_found++;
            points_found(number_of_points_found)=cell->node_index;distance_squared_of_points_found(number_of_points_found)=distance_squared_to_query_location;
            if(number_of_points_found==points_found.m) max_distance_squared=ARRAYS_COMPUTATIONS::Max(distance_squared_of_points_found);}
        else{
            if(number_of_points_found==points_found.m){ARRAYS_COMPUTATIONS::Heapify(distance_squared_of_points_found,points_found);number_of_points_found++;}
            // can't use heapify here, because we need to break at the position we place the new point
            int current_index=1;
            int left,right,index_of_largest;
            for(;;){
                left=2*current_index;right=2*current_index+1;index_of_largest=current_index;
                if(left>points_found.m) break;
                else if(right>points_found.m || distance_squared_of_points_found(left)>distance_squared_of_points_found(right)) index_of_largest=left;
                else index_of_largest=right;
                if(distance_squared_to_query_location>distance_squared_of_points_found(index_of_largest)) break; // we found a place to insert the new point
                exchange(distance_squared_of_points_found(current_index),distance_squared_of_points_found(index_of_largest));
                exchange(points_found(current_index),points_found(index_of_largest));
                current_index=index_of_largest;}
            distance_squared_of_points_found(current_index)=distance_squared_to_query_location;
            points_found(current_index)=cell->node_index;
            max_distance_squared=distance_squared_of_points_found(1);}} // copy the new bounding distance
}
//#####################################################################
template class KD_TREE<VECTOR<float,1> >;
template class KD_TREE<VECTOR<float,2> >;
template class KD_TREE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KD_TREE<VECTOR<double,1> >;
template class KD_TREE<VECTOR<double,2> >;
template class KD_TREE<VECTOR<double,3> >;
#endif
