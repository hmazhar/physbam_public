//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINTREE_CHILDREN  
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __BINTREE_CHILDREN__
#define __BINTREE_CHILDREN__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
namespace PhysBAM{ 

template<class T>
class BINTREE_CHILDREN:public NONCOPYABLE
{
public:
    BINTREE_CELL<T> children[2]; // the 2 cells stored in this structure
    BINTREE_CELL<T>* parent; // the cell for which this structure holds its 2 children 
    int* nodes; // the 3 nodes that the 2 children share - 0 if there are no nodes
    int* faces; // the 3 faces that the 2 children share - 0 if there is no face data
    VECTOR<T,1> parents_center; // the center of the cell that would contain these 2 children
    VECTOR<T,1> childrens_DX; // the size of one of the child cells
    int childrens_depth; // the depth of the children
private:
    static const int node_indices[2][2]; // indexed by [child_index][node_index]
    static const int face_indices[2][2]; // indexed by the [child_index][node_index]
public:

    BINTREE_CHILDREN(const bool use_nodes,const bool use_faces)
    {
        if(use_nodes) nodes=new int[3];else nodes=0; 
        if(use_faces) faces=new int[3];else faces=0;
    }

    ~BINTREE_CHILDREN()
    {delete[] nodes;delete[] faces;}

    int& Node(const int child_index,const int node_index)
    {assert(child_index >= 0 && child_index < 2);assert(node_index >= 0 && node_index < 2);assert(nodes);
    return nodes[node_indices[child_index][node_index]];};

    int Node(const int child_index,const int node_index) const
    {assert(child_index >= 0 && child_index < 2);assert(node_index >= 0 && node_index < 2);assert(nodes);
    return nodes[node_indices[child_index][node_index]];};

    int& Face(const int child_index,const int face_index)
    {assert(child_index >= 0 && child_index < 2);assert(face_index >= 0 && face_index < 2);assert(faces);
    return faces[face_indices[child_index][face_index]];};

    int Face(const int child_index,const int face_index) const
    {assert(child_index >= 0 && child_index < 2);assert(face_index >= 0 && face_index < 2);assert(faces);
    return faces[face_indices[child_index][face_index]];};

//#####################################################################
public:
    void Create_Face_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& face_count);
    void Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count);
    void Initialize_Root_Cell(int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,1>& center_root,const VECTOR<T,1>& DX);
    void Initialize_Pseudo_Root_Cells(int& number_of_cells,ARRAY<int,VECTOR<int,1> >& nodes_input,ARRAY<int,VECTOR<int,1> >& faces_input,const VECTOR<T,1>& center_input,const VECTOR<T,1>& childrens_DX_input);
    void Initialize_Non_Root_Cell(BINTREE_CELL<T>* parent_input,int& number_of_cells,ARRAY<BINTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_nodes,
        int& number_of_faces,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_faces,const VECTOR<T,1>& parents_center_in,const BINTREE_GRID<T>* grid);
//#####################################################################
};
}
#endif
#endif
