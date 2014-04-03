//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GREEN_CHILDREN_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __GREEN_CHILDREN_2D__
#define __GREEN_CHILDREN_2D__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
namespace PhysBAM{

class TRIANGLE_MESH;
template<class T> class RED_TRIANGLE;
template<class T> class RED_GREEN_GRID_2D;

template<class T>
class GREEN_CHILDREN_2D
{
public:
    RED_TRIANGLE<T>* parent; // the triangle for which this structure holds its children 
    ARRAY<VECTOR<int,3> > elements; // for each green child, its three nodes
    ARRAY<int> cells; // the index of each triangle
    int midpoints[3]; // the 3 midpoints nodes for the red parent, which can be 0 if they don't exist

    GREEN_CHILDREN_2D(RED_TRIANGLE<T>* parent_input=0)
        :parent(parent_input)
    {
        midpoints[0]=midpoints[1]=midpoints[2]=0;
    }

    ~GREEN_CHILDREN_2D()
    {}

    int& Midpoint(const int midpoint_index)
    {assert(midpoint_index >= 0 && midpoint_index < 3);
    return midpoints[midpoint_index];}
    
//#####################################################################
    void Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count);
    void Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count);
    bool Simulate_Green_Refine(int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid);
    void Finalize_Green_Refinement(int& number_of_cells);
    void Build_Triangle_Mesh(TRIANGLE_MESH& triangle_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_triangle_mapping,ARRAY<int>* node_to_particle_mapping) const;
    int Green_Leaf_Triangle(const VECTOR<T,2>& location) const;
    void Add_Ordered_Neighbors(ARRAY<ARRAY<int> >& neighbors,ARRAY<ARRAY<int> >& neighbor_links) const;
//#####################################################################
};
}
#endif
#endif
