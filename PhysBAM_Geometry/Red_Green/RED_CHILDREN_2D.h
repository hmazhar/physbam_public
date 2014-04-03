//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_CHILDREN_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __RED_CHILDREN_2D__
#define __RED_CHILDREN_2D__

#include <PhysBAM_Geometry/Red_Green/RED_TRIANGLE.h>
namespace PhysBAM{ 

template<class T>
class RED_CHILDREN_2D:public NONCOPYABLE
{
public:
    RED_TRIANGLE<T> children[4]; // the 4 triangles stored in this structure
    RED_TRIANGLE<T>* parent; // the triangle for which this structure holds its 4 children
    VECTOR<T,2> center; // the center of the parent triangle
    T childrens_size; // distance from center to Corner(2)
    int orientation; // orientation
    int nodes[6]; // the 6 nodes that the 4 children share
    int childrens_depth; // the depth of the children
private:
    static const int relative_orientations[4];
    static const int node_indices[4][3]; // indexed by [child_index][node_index]
public:

    RED_CHILDREN_2D()
    {}

    ~RED_CHILDREN_2D()
    {}
    
    void Clean_Memory()
    {for(int i=0;i<4;i++)children[i].Clean_Memory();}

    int& Node(const int child_index,const int node_index)
    {assert(child_index >= 0 && child_index < 4);assert(node_index >= 0 && node_index < 3);
    return nodes[node_indices[child_index][node_index]];};
    
    int& Corner(const int corner_index)
    {assert(corner_index >= 0 && corner_index < 3);
    return nodes[corner_index<<1];}
    
    int& Midpoint(const int midpoint_index)
    {assert(midpoint_index >= 0 && midpoint_index < 3);
    return nodes[(midpoint_index<<1)+1];}

    int Child_Orientation(const int child_index)
    {assert(child_index >= 0 && child_index < 4);
    return (orientation+relative_orientations[child_index])&3;}
    
//#####################################################################
    void Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count);
    void Initialize_Pseudo_Root_Triangle(int& number_of_cells,int node1,int node2,int node3,const VECTOR<T,2>& center_input,const T childrens_size_input,const int orientation_input);
    void Initialize_Non_Root_Triangle(RED_TRIANGLE<T>* parent_input,int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid);
//#####################################################################
};
}
#endif
#endif
