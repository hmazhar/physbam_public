//#####################################################################
// Copyright 2002-2006, Christopher Allocco, Robert Bridson, Kevin Der, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Duc Nguyen, Eftychios Sifakis, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_MESH
//##################################################################### 
#ifndef __TRIANGLE_MESH__
#define __TRIANGLE_MESH__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
namespace PhysBAM{

class SEGMENT_MESH;

class TRIANGLE_MESH:public SIMPLEX_MESH<2>
{
public:
    typedef VECTOR<int,3> ELEMENT_TYPE;
    typedef SIMPLEX_MESH<2> BASE;
private:
    typedef MESH_POLICY<BASE::dimension-1>::MESH T_BOUNDARY_MESH; // encourage this to eventually be "lifted" to SIMPLEX_MESH
public:

    using BASE::number_nodes;using BASE::elements;using BASE::neighbor_nodes;using BASE::incident_elements;using BASE::adjacent_elements;using BASE::neighbor_elements;
    using BASE::Simplex;using BASE::Node_In_Simplex;using BASE::Add_Connectivity;

    ARRAY<ARRAY<int> >* topologically_sorted_neighbor_nodes; // sorted version of neighbor_nodes for manifold meshes
    ARRAY<ARRAY<int> >* topologically_sorted_incident_elements; // sorted version of incident_elements for manifold meshes
    SEGMENT_MESH* segment_mesh; // segment mesh consisting of all the edges
    ARRAY<VECTOR<int,3> >* element_edges; // array of 3 indices for each triangle - edge element_edges(j,i) is edge j in triangle i
    ARRAY<ARRAY<int> >* edge_triangles; // for each edge, the indices of the incident triangles
    T_BOUNDARY_MESH* boundary_mesh;
    ARRAY<bool>* node_on_boundary; // whether or not a node is on the boundary
    ARRAY<int>* boundary_nodes;

    TRIANGLE_MESH(); // simplest constructor - null mesh
    TRIANGLE_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,3> >& triangle_list);
    TRIANGLE_MESH(const TRIANGLE_MESH& triangle_mesh);
    ~TRIANGLE_MESH();
    
    bool Node_In_Triangle(const int node,const int triangle_index) const
    {return Node_In_Simplex(node,triangle_index);}

    bool Segment_In_Triangle(const int segment_node_1,const int segment_node_2,const int triangle_index)
    {assert(segment_node_1!=segment_node_2);return Nodes_In_Simplex(VECTOR<int,2>(segment_node_1,segment_node_2),triangle_index);}

    void Replace_Node_In_Triangle(const int triangle,const int old_node,const int new_node)
    {Replace_Node_In_Simplex(triangle,old_node,new_node);}
    
    bool Edge_Neighbors(const int triangle1,const int triangle2) const
    {int i,j,k;elements(triangle1).Get(i,j,k);
    return Node_In_Triangle(i,triangle2)+Node_In_Triangle(j,triangle2)+Node_In_Triangle(k,triangle2) == 2;}

    int Other_Node(const int node1,const int node2,const int triangle) const
    {assert(Node_In_Triangle(node1,triangle) && Node_In_Triangle(node2,triangle));
    int i,j,k;elements(triangle).Get(i,j,k);return i^j^k^node1^node2;}
    
    void Other_Two_Nodes(const int node,const int triangle,int& other_node1,int& other_node2) const
    {assert(Node_In_Triangle(node,triangle));
    int i,j,k;elements(triangle).Get(i,j,k);
    if(i == node){other_node1=j;other_node2=k;}
    else if(j == node){other_node1=k;other_node2=i;}
    else{other_node1=i;other_node2=j;}} // k == node

    static bool Equivalent_Triangles(const int tri1_a,const int tri1_b,const int tri1_c,const int tri2_a,const int tri2_b,const int tri2_c)
    {if(tri1_a==tri2_a && ((tri1_b==tri2_b && tri1_c==tri2_c) || (tri1_b==tri2_c && tri1_c==tri2_b))) return true;
    if(tri1_a==tri2_b && ((tri1_b==tri2_c && tri1_c==tri2_a) || (tri1_b==tri2_a && tri1_c==tri2_c))) return true;
    if(tri1_a==tri2_c && ((tri1_b==tri2_a && tri1_c==tri2_b) || (tri1_b==tri2_b && tri1_c==tri2_a))) return true;
    return false;}

    static bool Equivalent_Oriented_Triangles(const int tri1_a,const int tri1_b,const int tri1_c,const int tri2_a,const int tri2_b,const int tri2_c)
    {return ((tri1_a==tri2_a && tri1_b==tri2_b && tri1_c==tri2_c) || (tri1_a==tri2_b && tri1_b==tri2_c && tri1_c==tri2_a) || (tri1_a==tri2_c && tri1_b==tri2_a && tri1_c==tri2_b));}

    static void Add_Ordered_Neighbors(ARRAY<int>& nodes,ARRAY<int>& links,const int neighbor1,const int neighbor2)
    {int index1,index2;
    if(!nodes.Find(neighbor1,index1)){nodes.Append(neighbor1);links.Append(0);index1=nodes.m;}
    if(!nodes.Find(neighbor2,index2)){nodes.Append(neighbor2);links.Append(0);index2=nodes.m;}
    links(index1)=index2;}

    static bool Face_Reversed_In_Simplex(const VECTOR<int,2>& segment_nodes,const VECTOR<int,3>& triangle_nodes)
    {int i=triangle_nodes.Find(segment_nodes[1]);if(triangle_nodes[i%3+1]==segment_nodes[2]) return false;
    assert(triangle_nodes[(i+1)%3+1]==segment_nodes[2]);return true;}

    SEGMENT_MESH& Get_Segment_Mesh()
    {if(!segment_mesh) Initialize_Segment_Mesh();return *segment_mesh;}

    int Triangle(const int node1,const int node2,const int node3) const
    {return Simplex(VECTOR<int,3>(node1,node2,node3));}

    void Mark_Edge_Connected_Component_Incident_On_A_Node(const int node,ARRAY<bool>& marked) const
    {Mark_Edge_Connected_Component_Incident_On_A_Node(node,1,marked);}

//#####################################################################
    void Delete_Auxiliary_Structures() PHYSBAM_OVERRIDE;
    void Refresh_Auxiliary_Structures() PHYSBAM_OVERRIDE;
    void Initialize_Square_Mesh(const int m,const int n,const bool reverse_triangles=false); // construct a regular m-by-n rectangular mesh
    void Initialize_Equilateral_Mesh(const int m,const int n); 
    void Initialize_Torus_Mesh(const int m,const int n); 
    void Initialize_Circle_Mesh(const int num_radial,const int num_tangential); // construct a circle
    void Initialize_Herring_Bone_Mesh(const int m,const int n); // construct a regular m-by-n herring bone rectangular mesh
    void Initialize_Cylinder_Mesh(const int m,const int n,const bool create_caps=true);
    void Initialize_Topologically_Sorted_Neighbor_Nodes();
    void Initialize_Topologically_Sorted_Incident_Elements();
    void Initialize_Segment_Mesh();
    void Initialize_Element_Edges();
    void Initialize_Edge_Triangles();
    void Initialize_Boundary_Mesh();
    void Initialize_Node_On_Boundary();
    void Initialize_Boundary_Nodes();
    void Non_Manifold_Nodes(ARRAY<int>& node_list);
private:
    void Mark_Edge_Connected_Component_Incident_On_A_Node(const int node,const int triangle_index_in_incident_triangles,ARRAY<bool>& marked) const;
public:
    int Adjacent_Triangle(const int triangle_index,const int node1,const int node2) const;
    int Triangles_On_Edge(const int node1,const int node2,ARRAY<int>* triangles_on_edge) const;
    bool Triangle_On_Edge(const int node1,const int node2) const;
    int Triangles_Across_Edge(const int triangle,const int node1,const int node2,ARRAY<int>& triangles_across_edge) const;
    void Make_Orientations_Consistent();
    bool Orientations_Consistent();
    void Identify_Connected_Components(ARRAY<int>& label);
    void Label_Connected_Component_With_ID(ARRAY<int>& label,const int triangle,const int id) const;
    bool Assert_Consistent() const PHYSBAM_OVERRIDE;
    void Set_Number_Nodes(const int number_nodes_input) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
