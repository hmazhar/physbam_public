//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLEX_MESH
//#####################################################################
#ifndef __SIMPLEX_MESH__
#define __SIMPLEX_MESH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class SEGMENT_MESH;

template<int d>
class SIMPLEX_MESH
{
public:
    enum WORKAROUND {dimension=d};

    int number_nodes; // number of nodes in the mesh
    ARRAY<VECTOR<int,d+1> > elements; // array of d+1 indices for each simplex - elements(i,t) is i'th index in simplex t
    ARRAY<ARRAY<int> >* neighbor_nodes; // for each node, list of neighboring nodes
    ARRAY<ARRAY<int> >* incident_elements; // for each node, list of neighboring simplices that contain it
    ARRAY<ARRAY<int> >* adjacent_elements; // for each simplex, list of (up to d+1) adjacent neighboring simplices
    ARRAY<ARRAY<int> >* neighbor_elements; // for each simplex, list of simplices sharing at least one node

    SIMPLEX_MESH(); // simplest constructor - null mesh
    SIMPLEX_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,d+1> >& simplex_list);
    SIMPLEX_MESH(const SIMPLEX_MESH& mesh);
    virtual ~SIMPLEX_MESH();

private:
    void operator=(const SIMPLEX_MESH&); // use Initialize_Mesh instead
public:

    void Initialize_Mesh(const int number_nodes_input,const ARRAY<VECTOR<int,d+1> >& simplex_list) // construct a mesh given a list of simplices
    {Clean_Memory();number_nodes=number_nodes_input;elements=simplex_list;}

    void Initialize_Mesh(const SIMPLEX_MESH& mesh) // works with the copy constructor
    {Initialize_Mesh(mesh.number_nodes,mesh.elements);}

    void Initialize_Mesh_With_Particle_Offset(const SIMPLEX_MESH& mesh,const int particle_offset)
    {elements.Resize(mesh.elements.m,false,false);
    for(int t=1;t<=elements.m;t++) elements(t)=mesh.elements(t)+particle_offset;
    number_nodes=0;} // TODO: currently leaves number_nodes uninitialized

    bool Node_In_Simplex(const int node,const int simplex) const
    {return elements(simplex).Contains(node);}

    template<int d2>
    bool Nodes_In_Simplex(const VECTOR<int,d2>& nodes,const int simplex) const
    {STATIC_ASSERT(d2<=d+1);const VECTOR<int,d+1>& element=elements(simplex);
    for(int i=1;i<=nodes.m;i++) if(!element.Contains(nodes[i])) return false;return true;}

    void Replace_Node_In_Simplex(const int simplex,const int old_node,const int new_node)
    {VECTOR<int,d+1>& element=elements(simplex);element[element.Find(old_node)]=new_node;}

    void Add_Nodes(const int new_nodes)
    {Set_Number_Nodes(number_nodes+new_nodes);}

    static int Other_Node(const VECTOR<int,d+1>& simplex_nodes,const VECTOR<int,d>& subsimplex_nodes)
    {for(int i=1;i<=d;i++) if(!subsimplex_nodes.Contains(simplex_nodes[i])) return simplex_nodes[i];
    assert(!subsimplex_nodes.Contains(simplex_nodes[d+1]));return simplex_nodes[d+1];}

    template<class T> static VECTOR<T,d> Node_Weights(const VECTOR<int,d+1>& simplex_nodes,const int node)
    {VECTOR<T,d> weights;int i=simplex_nodes.Find(node);assert(i);if(i<=d) weights(i)=(T)1;return weights;}

    template<class T,int d2> static VECTOR<VECTOR<T,d>,d2> Subsimplex_Weights(const VECTOR<int,d+1>& simplex_nodes,const VECTOR<int,d2>& subsimplex_nodes)
    {VECTOR<VECTOR<T,d>,d2> all_weights;for(int i=1;i<=d2;i++) all_weights(i)=Node_Weights<T>(simplex_nodes,subsimplex_nodes[i]);return all_weights;}

    template<class T_CONNECTIVITY> void Add_Connectivity(T_CONNECTIVITY& particle_connectivity) const
    {for(int t=1;t<=elements.m;t++) particle_connectivity.Union(elements(t));}

//#####################################################################
    virtual void Clean_Memory();
    virtual void Delete_Auxiliary_Structures();
    virtual void Refresh_Auxiliary_Structures();
    int Simplex(const VECTOR<int,d+1>& nodes) const;
    void Initialize_Neighbor_Nodes();
    void Initialize_Incident_Elements();
    void Initialize_Adjacent_Elements();
private: // helper function for Initialize_Adjacent_Elements
    void Find_And_Append_Adjacent_Elements(const int simplex,const VECTOR<int,d>& face);
public:
    void Initialize_Neighbor_Elements();
    int Add_Element_If_Not_Already_There(const VECTOR<int,d+1>& nodes);
    int Delete_Elements_With_Missing_Nodes(); // returns the number deleted
    void Delete_Sorted_Elements(const ARRAY<int>& deletion_list);
    void Delete_Sorted_Elements(const ARRAY<int>& deletion_list,HASHTABLE<int,int>& index_map);
    void Delete_Elements(ARRAY<int> deletion_list);
    int Number_Of_Nodes_With_Minimum_Valence();
    int Minimum_Valence(int* index=0); // ignores zero valence particles
    int Maximum_Valence(int* index=0);
    void Update_Adjacent_Elements_From_Incident_Elements(const int node);
    void Update_Neighbor_Nodes_From_Incident_Elements(const int node);
    virtual bool Assert_Consistent() const;
    virtual void Set_Number_Nodes(const int number_nodes_input);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const;
    template<class T> void Mark_Nodes_Referenced(ARRAY<T>& marks,const T& mark) const;
    template<int d2> void Simplices_On_Subsimplex(const VECTOR<int,d2>& subsimplex_nodes,ARRAY<int>& simplices_on_subsimplex) const;
//#####################################################################
};
}
#endif
