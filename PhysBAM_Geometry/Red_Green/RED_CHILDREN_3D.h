//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_CHILDREN_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __RED_CHILDREN_3D__
#define __RED_CHILDREN_3D__

#include <PhysBAM_Geometry/Red_Green/RED_TETRAHEDRON.h>
#include <PhysBAM_Geometry/Red_Green/TETRAHEDRAL_GROUP.h>
namespace PhysBAM{ 

template<class T>
class RED_CHILDREN_3D:public NONCOPYABLE
{
public:
    RED_TETRAHEDRON<T> children[8]; // the 8 tetrahedrons stored in this structure
    RED_TETRAHEDRON<T>* parent; // the tetrahedron for which this structure holds its 8 children
    VECTOR<T,3> center; // the center of the parent tetrahedron
    T childrens_size; // distance from center to midpoint of edge 01 or 23
    TETRAHEDRAL_GROUP<T> orientation; // orientation
    int nodes[10]; // the 10 nodes that the 4 children share
    int childrens_depth; // the depth of the children
private:
    static const TETRAHEDRAL_GROUP<T> relative_orientations[8];
    static const int node_indices[8][4]; // indexed by [child_index][node_index]
public:

    RED_CHILDREN_3D()
    {}

    ~RED_CHILDREN_3D()
    {}
    
    void Clean_Memory()
    {for(int i=0;i<8;i++)children[i].Clean_Memory();}

    int& Node(const int child_index,const int node_index)
    {assert(child_index >= 0 && child_index < 8);assert(node_index >= 0 && node_index < 4);
    return nodes[node_indices[child_index][node_index]];};
    
    int& Corner(const int corner_index)
    {assert(corner_index >= 0 && corner_index < 4);
    return nodes[corner_index];}
    
    int& Midpoint(const int midpoint_index)
    {assert(midpoint_index >= 0 && midpoint_index < 6);
    return nodes[midpoint_index+4];}

    TETRAHEDRAL_GROUP<T> Child_Orientation(const int child_index)
    {assert(child_index >= 0 && child_index < 8);
    return orientation*relative_orientations[child_index];}
    
//#####################################################################
    void Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count);
    void Initialize_Pseudo_Root_Tetrahedron(int& number_of_cells,int node1,int node2,int node3,int node4,const VECTOR<T,3>& center_input,const T childrens_size_input,const TETRAHEDRAL_GROUP<T>& orientation_input);
    void Initialize_Non_Root_Tetrahedron(RED_TETRAHEDRON<T>* parent_input,int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid);
//#####################################################################
};
}
#endif
#endif
