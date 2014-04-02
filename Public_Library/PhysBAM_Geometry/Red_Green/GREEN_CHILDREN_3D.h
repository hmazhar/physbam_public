//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GREEN_CHILDREN_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __GREEN_CHILDREN_3D__
#define __GREEN_CHILDREN_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
namespace PhysBAM{ 

class TETRAHEDRON_MESH;
template<class T> class RED_TETRAHEDRON;
template<class T> class RED_GREEN_GRID_3D;

template<class T>
class GREEN_CHILDREN_3D
{
public:
    RED_TETRAHEDRON<T>* parent; // the tetrahedron for which this structure holds its children 
    ARRAY<VECTOR<int,4> > elements; // for each green child, its four nodes
    ARRAY<int> cells; // the index of each tetrahedron
    int midpoints[6]; // the 6 midpoints nodes for the red parent, which can be 0 if they don't exist

    GREEN_CHILDREN_3D(RED_TETRAHEDRON<T>* parent_input=0)
        :parent(parent_input)
    {
        for(int i=0;i<6;i++)midpoints[i]=0;
    }

    ~GREEN_CHILDREN_3D()
    {}

    int& Midpoint(const int midpoint_index)
    {assert(midpoint_index >= 0 && midpoint_index < 6);
    return midpoints[midpoint_index];}
    
//#####################################################################
    void Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count);
    void Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count);
    bool Simulate_Green_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid);
    void Finalize_Green_Refinement(int& number_of_cells);
    void Build_Tetrahedron_Mesh(TETRAHEDRON_MESH& tetrahedron_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_tetrahedron_mapping,ARRAY<int>* node_to_particle_mapping) const;
    int Green_Leaf_Tetrahedron(const VECTOR<T,3>& location,const ARRAY<VECTOR<T,3> >& node_locations) const;
//#####################################################################
};
}
#endif
#endif
