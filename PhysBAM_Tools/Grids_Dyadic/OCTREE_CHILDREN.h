#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OCTREE_CHILDREN  
//#####################################################################
#ifndef __OCTREE_CHILDREN__
#define __OCTREE_CHILDREN__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
namespace PhysBAM{ 

template<class T>
class OCTREE_CHILDREN:public NONCOPYABLE
{
public:
    OCTREE_CELL<T> children[8]; // the 8 cells stored in this structure
    OCTREE_CELL<T>* parent; // the cell for which this structure holds its 8 children 
    int* nodes; // the 27 nodes that the 8 children share - 0 if there are no nodes
    int* faces; // the 36 faces that the 8 children share - 0 if there is no face data
    VECTOR<T,3> parents_center; // the center of the cell that would contain these 8 children
    VECTOR<T,3> childrens_DX; // the size of one of the child cells
    int childrens_depth; // the depth of the children
private:
    static const int node_indices[8][8]; // indexed by [child_index][node_index]
    static const int face_indices[8][6]; // indexed by the [child_index][node_index]
public:

    OCTREE_CHILDREN(const bool use_nodes,const bool use_faces)
    {
        if(use_nodes) nodes=new int[27];else nodes=0; 
        if(use_faces) faces=new int[36];else faces=0;
    }

    ~OCTREE_CHILDREN()
    {delete[] nodes;delete[] faces;}

    int& Node(const int child_index,const int node_index)
    {assert(child_index >= 0 && child_index < 8);assert(node_index >= 0 && node_index < 8);assert(nodes);
    return nodes[node_indices[child_index][node_index]];};

    int Node(const int child_index,const int node_index) const
    {assert(child_index >= 0 && child_index < 8);assert(node_index >= 0 && node_index < 8);assert(nodes);
    return nodes[node_indices[child_index][node_index]];};

    int& Face(const int child_index,const int face_index)
    {assert(child_index >= 0 && child_index < 8);assert(face_index >= 0 && face_index < 6);assert(faces);
    return faces[face_indices[child_index][face_index]];};

    int Face(const int child_index,const int face_index) const
    {assert(child_index >= 0 && child_index < 8);assert(face_index >= 0 && face_index < 6);assert(faces);
    return faces[face_indices[child_index][face_index]];};

private:
    static int Node_Index(const int i,const int j,const int k) // translates from 3d coordinates to the flat array that it is stored in
    {assert(i>=0 && i<3);assert(j>=0 && j<3);assert(k>=0 && k<3);return 9*k+3*j+i;}

//#####################################################################
public:
    void Create_Face_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& face_count);
    void Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count);
    void Initialize_Root_Cell(int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,3>& center_root,const VECTOR<T,3>& DX);
    void Initialize_Pseudo_Root_Cells(int& number_of_cells,ARRAY<int,VECTOR<int,1> >& nodes_input,ARRAY<int,VECTOR<int,1> >& faces_input,const VECTOR<T,3>& center_input,const VECTOR<T,3>& childrens_DX_input);
    void Initialize_Non_Root_Cell(OCTREE_CELL<T>* parent_input,int& number_of_cells,ARRAY<OCTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes,
        int& number_of_faces,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_faces,const VECTOR<T,3>& parents_center_in,const OCTREE_GRID<T>* grid);
    int Storage_Requirement()const;
private:
    void Edge_Neighbor_Check(int& number_of_nodes,const int x_offset,const int y_offset,const int z_offset,const int node_index,const int neighbor_node_index,const int child_number,
        const int node_number,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes,const OCTREE_GRID<T>* grid);
//#####################################################################
};
}
#endif
#endif
