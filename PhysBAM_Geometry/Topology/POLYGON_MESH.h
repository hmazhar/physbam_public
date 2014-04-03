//#####################################################################
// Copyright 2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGON_MESH
//#####################################################################
#ifndef __POLYGON_MESH__
#define __POLYGON_MESH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
namespace PhysBAM{

class SEGMENT_MESH;

class POLYGON_MESH
{
public:
    int number_nodes; // number of nodes in the mesh
    ARRAY<ARRAY<ARRAY<int> > > elements; // for each polygon, loop component, and vertex 
    SEGMENT_MESH* segment_mesh; // elements guaranteed to be lexographically ordered
    ARRAY<ARRAY<ARRAY<PAIR<int,bool> > > >* element_oriented_edges; // true indicates that the element has the edge oriented as in the segment mesh
    ARRAY<ARRAY<int> >* edge_elements; // for each element in segment mesh indicates which polygon elements contain it

    POLYGON_MESH(); // simplest constructor - null mesh
    POLYGON_MESH(const int number_nodes_input,const ARRAY<ARRAY<ARRAY<int> > >& polygon_list);
    POLYGON_MESH(const POLYGON_MESH& mesh);
    virtual ~POLYGON_MESH();

private:
    void operator=(const POLYGON_MESH&); // use Initialize_Mesh instead
public:

    void Initialize_Mesh(const int number_nodes_input,const ARRAY<ARRAY<ARRAY<int> > >& polygon_list) // construct a mesh given a list of polygons
    {Clean_Memory();number_nodes=number_nodes_input;elements=polygon_list;}
    
    void Initialize_Mesh(const POLYGON_MESH& mesh) // works with the copy constructor
    {Initialize_Mesh(mesh.number_nodes,mesh.elements);}
    
    void Add_Nodes(const int new_nodes)
    {Set_Number_Nodes(number_nodes+new_nodes);}

//#####################################################################
    virtual void Clean_Memory();
    virtual void Delete_Auxiliary_Structures();
    virtual void Refresh_Auxiliary_Structures();
    virtual bool Assert_Consistent() const;
    virtual void Set_Number_Nodes(const int number_nodes_input);
    void Initialize_Segment_Mesh();
    void Initialize_Element_Oriented_Edges();
    void Initialize_Edge_Elements();
    bool Oriented_Edge_In_Element(const int node1,const int node2,const int element) const;
    int Elements_On_Edge(const int node1,const int node2,ARRAY<int>* elements_on_edge) const;
    int Elements_On_Oriented_Edge(const int node1,const int node2,ARRAY<int>* elements_on_oriented_edge) const;
    int Opposite_Oriented_Element(const int element) const;
    void Split_Polygon_Edge(const int node1,const int node2,const int new_node);
//#####################################################################
};
}
#endif
