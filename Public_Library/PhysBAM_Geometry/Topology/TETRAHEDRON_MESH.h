//#####################################################################
// Copyright 2002-2006, Zhaosheng Bao, Robert Bridson, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Neil Molino, Avi Robinson-Mosher, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRON_MESH
//##################################################################### 
#ifndef __TETRAHEDRON_MESH__
#define __TETRAHEDRON_MESH__

#include <PhysBAM_Tools/Math_Tools/permutation.h>
#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
namespace PhysBAM{ 

class SEGMENT_MESH;
class TRIANGLE_MESH;

class TETRAHEDRON_MESH:public SIMPLEX_MESH<3>
{
public:
    typedef VECTOR<int,4> ELEMENT_TYPE;
    typedef SIMPLEX_MESH<3> BASE;
    using BASE::number_nodes;using BASE::elements;using BASE::neighbor_nodes;using BASE::incident_elements;using BASE::adjacent_elements;using BASE::neighbor_elements;
    using BASE::Simplex;using BASE::Node_In_Simplex;using BASE::Add_Connectivity;

    SEGMENT_MESH* segment_mesh; // segment mesh consisting of all the edges
    TRIANGLE_MESH* triangle_mesh; // triangle mesh consisting of all the triangles
    ARRAY<VECTOR<int,6> >* element_edges; // array of 6 indices for each tetrahedron - edges(j,i) is edge j in tetrahedron i
    ARRAY<VECTOR<int,4> >* tetrahedron_faces; // array of 4 indices for each tetrahedron - faces(j,i) is the face opposite vertex j in tetrahedron i
    TRIANGLE_MESH* boundary_mesh;
    ARRAY<bool>* node_on_boundary; // node_on_boundary(i) is true if node i is on the boundary
    ARRAY<int>* boundary_nodes; // array containing indices of all particles on the boundary  
    ARRAY<ARRAY<int> >* edge_tetrahedrons; // for each edge, list of tets containing that edge
    ARRAY<VECTOR<int,2> >* triangle_tetrahedrons; // for each face, list of incident tets
private:
    bool owns_segment_mesh;
public:

    TETRAHEDRON_MESH(); // simplest constructor - null mesh
    TETRAHEDRON_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,4> >& tetrahedron_list);
    TETRAHEDRON_MESH(const TETRAHEDRON_MESH& tetrahedron_mesh);
    ~TETRAHEDRON_MESH();

    bool Node_In_Tetrahedron(const int node,const int tetrahedron) const
    {return Node_In_Simplex(node,tetrahedron);}
    
    bool Triangle_In_Tetrahedron(const int node1,const int node2,const int node3,const int tetrahedron) const
    {return Nodes_In_Simplex(VECTOR<int,3>(node1,node2,node3),tetrahedron);}

    bool Face_Neighbors(const int tetrahedron1,const int tetrahedron2) const
    {int i,j,k,l;elements(tetrahedron1).Get(i,j,k,l);
    return Triangle_In_Tetrahedron(j,k,l,tetrahedron2) || Triangle_In_Tetrahedron(i,k,l,tetrahedron2) || 
        Triangle_In_Tetrahedron(i,j,l,tetrahedron2) || Triangle_In_Tetrahedron(i,j,k,tetrahedron2);}

    void Other_Three_Nodes(const int node,const int tetrahedron,int& other_node1,int& other_node2,int& other_node3) const
    {assert(Node_In_Tetrahedron(node,tetrahedron));
    int i,j,k,l;elements(tetrahedron).Get(i,j,k,l);
    if(i == node){other_node1=j;other_node2=k;other_node3=l;}
    else if(j == node){other_node1=i;other_node2=l;other_node3=k;}
    else if(k == node){other_node1=l;other_node2=i;other_node3=j;}
    else{other_node1=k;other_node2=j;other_node3=i;}} // l == node

    void Other_Two_Nodes(const int node1,const int node2,const int tetrahedron,int& other_node1,int& other_node2) const
    {assert(Node_In_Tetrahedron(node1,tetrahedron) && Node_In_Tetrahedron(node2,tetrahedron));
    int i,j,k,l;elements(tetrahedron).Get(i,j,k,l);
    if(i == node1){if(j == node2){other_node1=k;other_node2=l;}else if(k == node2){other_node1=l;other_node2=j;}else{other_node1=j;other_node2=k;}}
    else if(j == node1){if(i == node2){other_node1=l;other_node2=k;}else if(k == node2){other_node1=i;other_node2=l;}else{other_node1=k;other_node2=i;}}
    else if(k == node1){if(i == node2){other_node1=j;other_node2=l;}else if(l == node2){other_node1=i;other_node2=j;}else{other_node1=l;other_node2=i;}}
    else{if(i == node2){other_node1=k;other_node2=j;}else if(k == node2){other_node1=j;other_node2=i;}else{other_node1=i;other_node2=k;}}} // l == node1

    SEGMENT_MESH& Get_Segment_Mesh()
    {if(!segment_mesh) Initialize_Segment_Mesh();return *segment_mesh;}

    int Tetrahedron(const int node1,const int node2,const int node3,const int node4) const
    {return Simplex(VECTOR<int,4>(node1,node2,node3,node4));}

    void Mark_Face_Connected_Component_Incident_On_A_Node(const int node,ARRAY<bool>& marked) const
    {Mark_Face_Connected_Component_Incident_On_A_Node(node,1,marked);}

    void Tetrahedrons_On_Edge(const VECTOR<int,2>& nodes,ARRAY<int>& tetrahedrons_on_edge) const
    {Simplices_On_Subsimplex(nodes,tetrahedrons_on_edge);}

    void Tetrahedrons_On_Face(const VECTOR<int,3>& nodes,ARRAY<int>& tetrahedrons_on_face) const
    {Simplices_On_Subsimplex(nodes,tetrahedrons_on_face);}

    static bool Face_Reversed_In_Simplex(const VECTOR<int,3>& face_nodes,const VECTOR<int,4>& tet_nodes)
    {for(int i=1;i<=24;i++) if(permute_four(tet_nodes,i).Remove_Index(4)==face_nodes) return !permutation_of_four_is_even(i);
    PHYSBAM_FATAL_ERROR();}

//#####################################################################
    void Delete_Auxiliary_Structures() PHYSBAM_OVERRIDE;
    void Refresh_Auxiliary_Structures() PHYSBAM_OVERRIDE;
    void Initialize_Octahedron_Mesh(const int m,const int n,const int p);
    void Initialize_Cube_Mesh(const int m,const int n,const int p); // 5 tetrahedrons per cube
    void Initialize_Prismatic_Cube_Mesh(const int m,const int n,const int p); // 6 tetrahedra per cube
private:
    int Number_Of_Tetrahedrons_Across_Face(const int tetrahedron,const int node1,const int node2,const int node3) const;
public:
    void Initialize_Segment_Mesh();
    void Initialize_Segment_Mesh_From_Triangle_Mesh();
    void Initialize_Triangle_Mesh();
    void Initialize_Element_Edges();
    void Initialize_Tetrahedron_Faces();
    void Initialize_Boundary_Mesh();
    void Initialize_Node_On_Boundary();
    void Initialize_Boundary_Nodes();
    void Initialize_Edge_Tetrahedrons();
    void Initialize_Triangle_Tetrahedrons();
    int Number_Of_Boundary_Tetrahedrons();
    int Number_Of_Interior_Tetrahedrons();
    bool Edge_Neighbors(const int tet1,const int tet2) const;
    int Number_Of_Edge_Neighbors(const int segment) const;
    void Initialize_Segment_Mesh_Of_Subset(SEGMENT_MESH& segment_mesh_of_subset,const ARRAY<bool>& subset) const;
    void Initialize_Boundary_Mesh_Of_Subset(TRIANGLE_MESH& boundary_mesh_of_subset,const ARRAY<bool>& subset);
    void Initialize_Boundary_Mesh_With_T_Junctions(TRIANGLE_MESH& boundary_mesh_with_t_junctions,const ARRAY<int>& t_junctions,const ARRAY<VECTOR<int,2> >& t_junction_parents);
    void Initialize_Bending_Tetrahedrons(TRIANGLE_MESH& triangle_mesh);
private:
    void Mark_Face_Connected_Component_Incident_On_A_Node(const int node,const int tetrahedron_index_in_incident_tetrahedrons,ARRAY<bool>& marked) const;
public:
    int Tetrahedrons_Across_Face(const int tetrahedron,const int node1,const int node2,const int node3,ARRAY<int>& tetrahedrons_across_face) const;
    void Identify_Face_Connected_Components(ARRAY<int>& label);
    void Identify_Edge_Connected_Components(ARRAY<int>& label);
    bool Assert_Consistent() const PHYSBAM_OVERRIDE;
    void Set_Number_Nodes(const int number_nodes_input) PHYSBAM_OVERRIDE;
//#####################################################################
};   
}
#endif
