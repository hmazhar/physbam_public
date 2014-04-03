//#####################################################################
// Copyright 2004-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KD_TREE
//#####################################################################
#ifndef __KD_TREE__
#define __KD_TREE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE_NODE.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/POINTER_POOL.h>
namespace PhysBAM{

template<class TV>
class KD_TREE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    KD_TREE_NODE<T>* root_node;
    bool store_values_on_internal_nodes;
    POINTER_POOL<KD_TREE_NODE<T> > pool;

    KD_TREE(const bool store_values_on_internal_nodes_input=true)
        :root_node(0),store_values_on_internal_nodes(store_values_on_internal_nodes_input)
    {}

    ~KD_TREE()
    {}

private:
    int Choose_Partition_Index_Using_Internal_Nodes(const int first_index,const int last_index) const 
    {int elements=last_index-first_index+1,filled_subtree_elements=(1<<(integer_log(elements+1)-1))-1;
    return first_index+filled_subtree_elements+min(elements-2*filled_subtree_elements-1,filled_subtree_elements+1);} // this index goes to the root of the subtree

    int Choose_Partition_Index_Using_Leaf_Nodes(const int first_index,const int last_index) const
    {int elements=last_index-first_index+1,filled_subtree_elements=(1<<integer_log(elements));
    return first_index+min(elements-filled_subtree_elements/2,filled_subtree_elements);} // this index goes to the right subtree

    int Choose_Partition_Axis(const TV& DX) const
    {return DX.Arg_Max();}

    class Partition_Helper_Less{
    public:
        ARRAY_VIEW<const TV>* points;int axis;T split_value;
        Partition_Helper_Less(ARRAY_VIEW<const TV>* p,int a,T s){points=p;axis=a;split_value=s;}
        bool operator()(int i){return (*points)(i)[axis]<split_value;}};

    class Partition_Helper_Less_Equal{
    public:
        ARRAY_VIEW<const TV>* points;int axis;T split_value;
        Partition_Helper_Less_Equal(ARRAY_VIEW<const TV>* p,int a,T s){points=p;axis=a;split_value=s;}
        bool operator()(int i){return (*points)(i)[axis]<=split_value;}};

public:

//#####################################################################
    KD_TREE_NODE<T>* Leaf_Node(const TV& point) const;
    void Create_KD_Tree(ARRAY_VIEW<const TV> nodes);
    void Create_Left_Balanced_KD_Tree(ARRAY_VIEW<const TV> nodes_to_balance);
    void Create_Left_Balanced_KD_Tree_With_Grouping(ARRAY_VIEW<const TV> nodes_to_balance,ARRAY<ARRAY<int> >& points_in_group,const int max_points_in_group);
    void Find_Points_Within_Radius(const TV& location,const T max_distance_squared,ARRAY<int>& points_found,ARRAY<T>& distance_squared_of_points_found,
        ARRAY_VIEW<const TV> original_points) const;
    void Locate_Nearest_Neighbors(const TV& location,const T max_distance_squared,ARRAY<int>& points_found,ARRAY<T>& distance_squared_of_points_found,
        int& number_points_found,T& max_distance_squared_of_found_points,ARRAY_VIEW<const TV> original_points) const;
private:
    void Find_Points_Within_Radius_Helper(const KD_TREE_NODE<T>* cell,const TV& location,const T max_distance_squared,ARRAY<int>& points_found,
        ARRAY<T>& distance_squared_of_points_found,ARRAY_VIEW<const TV> original_points) const;
    void Locate_Nearest_Neighbors_Helper(const KD_TREE_NODE<T>* cell,const TV& location,T& max_distance_squared,int& number_of_points_found,ARRAY<int>& points_found,
        ARRAY<T>& distance_squared_of_points_found,ARRAY_VIEW<const TV> original_points) const;
    KD_TREE_NODE<T>* Leaf_Node_Helper(const TV& point,KD_TREE_NODE<T>* node) const;
    void Balance_Sub_KD_Tree_Using_Internal_Nodes(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,
        RANGE<TV>& box);
    void Balance_Sub_KD_Tree_Using_Leaf_Nodes(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,
        RANGE<TV>& box);
    void Balance_Sub_KD_Tree_Using_Leaf_Nodes_With_Grouping(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,
        RANGE<TV>& box,ARRAY<int>& point_group,int& current_group_index,const int max_points_in_group);
    void Sub_KD_Tree_Using_Leaf_Nodes(KD_TREE_NODE<T>* cell,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,RANGE<TV>& box);
    void Median_Split(const int required_partition_index,const int first_index,const int last_index,ARRAY_VIEW<const TV> points,ARRAY_VIEW<int> permutation_array,const int axis);
//#####################################################################
};
}
#endif
