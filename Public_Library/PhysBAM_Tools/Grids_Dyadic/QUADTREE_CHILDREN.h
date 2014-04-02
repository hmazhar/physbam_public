//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADTREE_CHILDREN  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __QUADTREE_CHILDREN__
#define __QUADTREE_CHILDREN__

#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
namespace PhysBAM{ 

template<class T>
class QUADTREE_CHILDREN:public NONCOPYABLE
{
public:
    QUADTREE_CELL<T> children[4]; // the 4 cells stored in this structure
    QUADTREE_CELL<T>* parent; // the cell for which this structure holds its 4 children 
    int* nodes; // the 9 nodes that the 4 children share - 0 if there are no nodes
    int* faces; // the 12 faces that the 4 children share - 0 if there is no face data
    VECTOR<T,2> parents_center; // the center of the cell that would contain these 4 children
    VECTOR<T,2> childrens_DX; // the size of one of the child cells
    int childrens_depth; // the depth of the children
private:
    static const int node_indices[4][4]; // indexed by [child_index][node_index]
    static const int face_indices[4][4]; // indexed by the [child_index][node_index]
public:

    QUADTREE_CHILDREN(const bool use_nodes,const bool use_faces)
    {
        if(use_nodes) nodes=new int[9];else nodes=0; 
        if(use_faces) faces=new int[12];else faces=0;
    }

    ~QUADTREE_CHILDREN()
    {delete[] nodes;delete[] faces;}

    int& Node(const int child_index,const int node_index)
    {assert(child_index >= 0 && child_index < 4);assert(node_index >= 0 && node_index < 4);assert(nodes);
    return nodes[node_indices[child_index][node_index]];};

    int Node(const int child_index,const int node_index) const
    {assert(child_index >= 0 && child_index < 4);assert(node_index >= 0 && node_index < 4);assert(nodes);
    return nodes[node_indices[child_index][node_index]];};

    int& Face(const int child_index,const int face_index)
    {assert(child_index >= 0 && child_index < 4);assert(face_index >= 0 && face_index < 4);assert(faces);
    return faces[face_indices[child_index][face_index]];};

    int Face(const int child_index,const int face_index) const
    {assert(child_index >= 0 && child_index < 4);assert(face_index >= 0 && face_index < 4);assert(faces);
    return faces[face_indices[child_index][face_index]];};

private:
    static int Node_Index(const int i,const int j) // translates from 2d coordinates to the flat array that it is stored in
    {assert(i>=0 && i<3);assert(j>=0 && j<3);return 3*j+i;}

//#####################################################################
public:
    void Create_Face_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& face_count);
    void Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count);
    void Initialize_Root_Cell(int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,2>& center_root,const VECTOR<T,2>& DX);
    void Initialize_Pseudo_Root_Cells(int& number_of_cells,ARRAY<int,VECTOR<int,1> >& nodes_input,ARRAY<int,VECTOR<int,1> >& faces_input,const VECTOR<T,2>& center_input,const VECTOR<T,2>& childrens_DX_input);
    void Initialize_Non_Root_Cell(QUADTREE_CELL<T>* parent_input,int& number_of_cells,ARRAY<QUADTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_nodes,
        int& number_of_faces,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_faces,const VECTOR<T,2>& parents_center_in,const QUADTREE_GRID<T>* grid);
//#####################################################################
};
}
#endif
#endif
