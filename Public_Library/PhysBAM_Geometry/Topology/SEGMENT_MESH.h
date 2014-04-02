//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Neil Molino, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_MESH
//#####################################################################
#ifndef __SEGMENT_MESH__
#define __SEGMENT_MESH__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
namespace PhysBAM{

class SEGMENT_MESH:public SIMPLEX_MESH<1>
{
public:
    typedef VECTOR<int,2> ELEMENT_TYPE;
    typedef SIMPLEX_MESH<1> BASE;
    typedef MESH_POLICY<BASE::dimension-1>::MESH T_BOUNDARY_MESH; // encourage this to eventually be "lifted" to SIMPLEX_MESH
    using BASE::number_nodes;using BASE::elements;using BASE::neighbor_nodes;using BASE::incident_elements;
    using BASE::Simplex;using BASE::Node_In_Simplex;using BASE::Replace_Node_In_Simplex;using BASE::Add_Connectivity;

    ARRAY<ARRAY<VECTOR<int,2> > >* connected_segments;  // a list of connected segment components, in no particular order
    ARRAY<ARRAY<int> >* ordered_loop_nodes;  // an ordered list of adjacent nodes for each connected component that is a simple loop (useful for boundaries)
    
    T_BOUNDARY_MESH* boundary_mesh;

    SEGMENT_MESH();
    SEGMENT_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,2> >& segment_list);
    SEGMENT_MESH(const SEGMENT_MESH& segment_mesh);
    ~SEGMENT_MESH();
    
    int Get_Opposite_Endpoint(const int segment,const int node)
    {assert(Node_In_Segment(node,segment));if(elements(segment)(1) == node) return elements(segment)(2);else return elements(segment)(1);}

    int Segment(const int node1,const int node2) const
    {return Simplex(VECTOR<int,2>(node1,node2));}

    bool Node_In_Segment(const int node,const int segment) const
    {return Node_In_Simplex(node,segment);}

    bool Segments_Adjacent(const int segment1,const int segment2) const
    {int i,j;elements(segment1).Get(i,j);return Node_In_Segment(i,segment2) || Node_In_Segment(j,segment2);}

    void Replace_Node_In_Segment(const int segment,const int old_node,const int new_node)
    {Replace_Node_In_Simplex(segment,old_node,new_node);}

    SEGMENT_MESH& Get_Segment_Mesh()
    {return *this;}

//#####################################################################
    void Delete_Auxiliary_Structures() PHYSBAM_OVERRIDE;
    void Refresh_Auxiliary_Structures() PHYSBAM_OVERRIDE;
    void Initialize_Connected_Segments();
    void Initialize_Ordered_Loop_Nodes();
    void Initialize_Straight_Mesh(const int number_of_points,bool loop=false);
    void Initialize_Boundary_Mesh();
    bool Assert_Consistent() const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
