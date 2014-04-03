//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Eran Guendelman,
// Geoffrey Irving, Frank Losasso, Neil Molino, Duc Nguyen, Eftychios Sifakis, Jerry Talton, Joseph Teran, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/exchange_sort.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
TRIANGLE_MESH::
TRIANGLE_MESH() // simplest constructor - null mesh
    :topologically_sorted_neighbor_nodes(0),topologically_sorted_incident_elements(0),
    segment_mesh(0),element_edges(0),edge_triangles(0),boundary_mesh(0),node_on_boundary(0),boundary_nodes(0)
{}
//#####################################################################
// Constructor
//#####################################################################
TRIANGLE_MESH::
TRIANGLE_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,3> >& triangle_list)
    :SIMPLEX_MESH<2>(number_nodes_input,triangle_list),topologically_sorted_neighbor_nodes(0),topologically_sorted_incident_elements(0),
    segment_mesh(0),element_edges(0),edge_triangles(0),boundary_mesh(0),node_on_boundary(0),boundary_nodes(0)
{
    Initialize_Mesh(number_nodes_input,triangle_list);
}
//#####################################################################
// Constructor
//#####################################################################
TRIANGLE_MESH::
TRIANGLE_MESH(const TRIANGLE_MESH& triangle_mesh)
    :SIMPLEX_MESH<2>(triangle_mesh),topologically_sorted_neighbor_nodes(0),topologically_sorted_incident_elements(0),
    segment_mesh(0),element_edges(0),edge_triangles(0),boundary_mesh(0),node_on_boundary(0),boundary_nodes(0)
{
    Initialize_Mesh(triangle_mesh);
}
//#####################################################################
// Destructor
//#####################################################################
TRIANGLE_MESH::
~TRIANGLE_MESH()
{
    delete topologically_sorted_neighbor_nodes;delete topologically_sorted_incident_elements;delete segment_mesh;
    delete element_edges;delete edge_triangles;delete boundary_mesh;delete node_on_boundary;delete boundary_nodes;
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
void TRIANGLE_MESH::
Delete_Auxiliary_Structures()
{
    SIMPLEX_MESH<2>::Delete_Auxiliary_Structures();
    delete topologically_sorted_neighbor_nodes;topologically_sorted_neighbor_nodes=0;delete topologically_sorted_incident_elements;topologically_sorted_incident_elements=0;
    delete segment_mesh;segment_mesh=0;delete element_edges;element_edges=0;
    delete edge_triangles;edge_triangles=0;delete boundary_mesh;boundary_mesh=0;
    delete node_on_boundary;node_on_boundary=0;delete boundary_nodes;boundary_nodes=0;
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
void TRIANGLE_MESH::
Refresh_Auxiliary_Structures()
{
    SIMPLEX_MESH<2>::Refresh_Auxiliary_Structures();
    if(topologically_sorted_neighbor_nodes) Initialize_Topologically_Sorted_Neighbor_Nodes();
    if(topologically_sorted_incident_elements) Initialize_Topologically_Sorted_Incident_Elements();
    if(segment_mesh) Initialize_Segment_Mesh();if(element_edges) Initialize_Element_Edges();if(edge_triangles) Initialize_Edge_Triangles();
    if(boundary_mesh) Initialize_Boundary_Mesh();if(node_on_boundary) Initialize_Node_On_Boundary();if(boundary_nodes) Initialize_Boundary_Nodes();
}
//#####################################################################
// Function Initialize_Square_Mesh
//#####################################################################
void TRIANGLE_MESH::
Initialize_Square_Mesh(const int m,const int n,const bool reverse_triangles) // construct a regular m-by-n rectangular mesh
{
    Clean_Memory();number_nodes=m*n;elements.Exact_Resize(2*(m-1)*(n-1));
    if(reverse_triangles){
        int t=0;for(int i=1;i<=m-1;i++)for(int j=1;j<=n-1;j++){ // counterclockwise node ordering
            elements(++t).Set(i+m*(j-1),i+1+m*(j-1),i+m*j);elements(++t).Set(i+1+m*(j-1),i+1+m*j,i+m*j);}}
    else{
        int t=0;for(int i=1;i<=m-1;i++)for(int j=1;j<=n-1;j++){ // counterclockwise node ordering
            elements(++t).Set(i+m*(j-1),i+1+m*(j-1),i+1+m*j);elements(++t).Set(i+m*(j-1),i+1+m*j,i+m*j);}}
}
//#####################################################################
// Function Initialize_Equilateral_Mesh
//#####################################################################
void TRIANGLE_MESH::
Initialize_Equilateral_Mesh(const int m,const int n)
{
    Clean_Memory();number_nodes=m*n;elements.Exact_Resize(2*(m-1)*(n-1));
    int t=0;for(int i=1;i<=m-1;i++)for(int j=1;j<=n-1;j++){ // counterclockwise node ordering
        if(j%2){elements(++t).Set(i+m*(j-1),i+1+m*(j-1),i+m*j);elements(++t).Set(i+1+m*(j-1),i+1+m*j,i+m*j);}
        else{elements(++t).Set(i+m*(j-1),i+1+m*(j-1),i+1+m*j);elements(++t).Set(i+m*(j-1),i+1+m*j,i+m*j);}}
}
//#####################################################################
// Function Initialize_Torus_Mesh
//#####################################################################
void TRIANGLE_MESH::
Initialize_Torus_Mesh(const int m,const int n)
{
    Clean_Memory();number_nodes=m*n;elements.Exact_Resize(2*m*n);assert((n&1)==0);
    int t=0;for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){ // counterclockwise node ordering
        int i1=i==m?1:i+1,j1=j==n?1:j+1;
        if(j&1){elements(++t).Set(i+m*(j-1),i1+m*(j-1),i+m*(j1-1));elements(++t).Set(i1+m*(j-1),i1+m*(j1-1),i+m*(j1-1));}
        else{elements(++t).Set(i+m*(j-1),i1+m*(j-1),i1+m*(j1-1));elements(++t).Set(i+m*(j-1),i1+m*(j1-1),i+m*(j1-1));}}
}
//#####################################################################
// Function Initialize_Circle_Mesh
//#####################################################################
void TRIANGLE_MESH::
Initialize_Circle_Mesh(const int num_radial,const int num_tangential) // construct a circle
{
    Clean_Memory();number_nodes=num_radial*num_tangential;elements.Exact_Resize(2*(num_radial-1)*num_tangential);
    int t=0,n=num_tangential;for(int i=0;i<num_radial-1;i++)for(int j=1;j<=n;j++){ // counterclockwise node ordering
        if(j&1){elements(++t).Set(i*n+j,i*n+(j%n)+1,(i+1)*n+j);elements(++t).Set((i+1)*n+j,i*n+(j%n)+1,(i+1)*n+(j%n)+1);}
        else{elements(++t).Set(i*n+j,i*n+(j%n)+1,(i+1)*n+j%n+1);elements(++t).Set((i+1)*n+j,i*n+j,(i+1)*n+(j%n)+1);}}
}
//#####################################################################
// Function Initialize_Herring_Bone_Mesh
//#####################################################################
void TRIANGLE_MESH::
Initialize_Herring_Bone_Mesh(const int m,const int n) // construct a regular m-by-n herring bone rectangular mesh
{
    Clean_Memory();number_nodes=m*n;elements.Exact_Resize(2*(m-1)*(n-1));
    int t=0;for(int i=1;i<=m-1;i++)for(int j=1;j<=n-1;j++){ // counterclockwise node ordering
        if(i%2){elements(++t).Set(i+m*(j-1),i+1+m*(j-1),i+m*j);elements(++t).Set(i+1+m*(j-1),i+1+m*j,i+m*j);}
        else{elements(++t).Set(i+m*(j-1),i+1+m*(j-1),i+1+m*j);elements(++t).Set(i+m*(j-1),i+1+m*j,i+m*j);}}
}
//#####################################################################
// Function Initialize_Cylinder_Mesh
//#####################################################################
void TRIANGLE_MESH::
Initialize_Cylinder_Mesh(const int m,const int n,const bool create_caps)
{
    Clean_Memory();int t=0;
    if(create_caps){elements.Exact_Resize(2*m*n);number_nodes=m*n+2;}
    else{elements.Exact_Resize(2*(m-1)*n);number_nodes=m*n;}
    for(int j=1;j<=n;j++){int j_1=j==n?1:j+1;
        for(int i=1;i<=m-1;i++){elements(++t).Set(j+(i-1)*n,j+i*n,j_1+i*n);elements(++t).Set(j+(i-1)*n,j_1+i*n,j_1+(i-1)*n);}
        if(create_caps){elements(++t).Set(m*n+1,j,j_1);elements(++t).Set(m*n+2,j_1+(m-1)*n,j+(m-1)*n);}}
}
//#####################################################################
// Function Initialize_Topologically_Sorted_Neighbor_Nodes
//#####################################################################
// for each node, sort its list of neighbors in the same order that triangle corners are listed 
// we assume triangle's corners are listed in a topologically consistent fashion, e.g. all counterclockwise
void TRIANGLE_MESH::
Initialize_Topologically_Sorted_Neighbor_Nodes()
{
    delete topologically_sorted_neighbor_nodes;
    ARRAY<ARRAY<int> > neighbors(number_nodes),neighbor_links(number_nodes);
    for(int t=1;t<=elements.m;t++){
        int i,j,k;elements(t).Get(i,j,k);
        Add_Ordered_Neighbors(neighbors(i),neighbor_links(i),j,k);
        Add_Ordered_Neighbors(neighbors(j),neighbor_links(j),k,i);
        Add_Ordered_Neighbors(neighbors(k),neighbor_links(k),i,j);}
    topologically_sorted_neighbor_nodes=new ARRAY<ARRAY<int> >(number_nodes);
    for(int i=1;i<=number_nodes;i++) if(neighbors(i).m){
        (*topologically_sorted_neighbor_nodes)(i).Exact_Resize(neighbors(i).m);
        ARRAY<bool> not_first(neighbors(i).m);for(int j=1;j<=neighbors(i).m;j++) if(neighbor_links(i)(j)) not_first(neighbor_links(i)(j))=true;
        int node_index=1;while(node_index <= neighbors(i).m && not_first(node_index)) node_index++; // now find the first node in the linked list
        if(node_index > neighbors(i).m) node_index=1; // if we have a cycle (i is in the interior), just use 1
        for(int j=1;j<=neighbors(i).m;j++){(*topologically_sorted_neighbor_nodes)(i)(j)=neighbors(i)(node_index);node_index=neighbor_links(i)(node_index);}}
}
//#####################################################################
// Function Initialize_Topologically_Sorted_Incident_Elements
//#####################################################################
// we assume triangles corners are listed in a topologically consistent fashion, e.g. all counterclockwise
void TRIANGLE_MESH::
Initialize_Topologically_Sorted_Incident_Elements()
{
    delete topologically_sorted_incident_elements;topologically_sorted_incident_elements=new ARRAY<ARRAY<int> >(number_nodes);
    Initialize_Incident_Elements();
    bool topologically_sorted_neighbor_nodes_defined=topologically_sorted_neighbor_nodes!=0;if(!topologically_sorted_neighbor_nodes_defined) Initialize_Topologically_Sorted_Neighbor_Nodes();
    for(int p=1;p<=number_nodes;p++){
        ARRAY<int>& neighbors=(*topologically_sorted_neighbor_nodes)(p);
        int m=neighbors.m;assert(m>=2);
        int last_triangle=Triangle(p,neighbors(1),neighbors(m));
        if(last_triangle){(*topologically_sorted_incident_elements)(p).Exact_Resize(m);(*topologically_sorted_incident_elements)(p)(m)=last_triangle;}
        else (*topologically_sorted_incident_elements)(p).Exact_Resize(m-1);
        for(int a=1;a<m;a++){int t=Triangle(p,neighbors(a),neighbors(a+1));assert(t);(*topologically_sorted_incident_elements)(p)(a)=t;}}
    if(!topologically_sorted_neighbor_nodes_defined){delete topologically_sorted_neighbor_nodes;topologically_sorted_neighbor_nodes=0;}
}
//#####################################################################
// Function Initialize_Segment_Mesh
//#####################################################################
void TRIANGLE_MESH::
Initialize_Segment_Mesh()
{
    delete segment_mesh;
    bool neighbor_nodes_defined=neighbor_nodes!=0;if(!neighbor_nodes_defined) Initialize_Neighbor_Nodes();
    // number of edges = half the sum of the degree (or a little more if there are loop edges)
    int total_degree=0;for(int i=1;i<=number_nodes;i++) total_degree+=(*neighbor_nodes)(i).m;
    segment_mesh=new SEGMENT_MESH();
    segment_mesh->number_nodes=number_nodes;
    segment_mesh->elements.Preallocate((total_degree+1)/2); // add one for subtle optimization purposes
    for(int node=1;node<=neighbor_nodes->m;node++) for(int k=1;k<=(*neighbor_nodes)(node).m;k++) // do nodes in ascending order so that no edges are counted more than once
        if(node<=(*neighbor_nodes)(node)(k)) segment_mesh->elements.Append(VECTOR<int,2>(node,(*neighbor_nodes)(node)(k)));
    if(!neighbor_nodes_defined){delete neighbor_nodes;neighbor_nodes=0;}
}
//#####################################################################
// Function Initialize_Element_Edges
//#####################################################################
void TRIANGLE_MESH::
Initialize_Element_Edges()
{
    delete element_edges;element_edges=new ARRAY<VECTOR<int,3> >(elements.m);
    if(!segment_mesh) Initialize_Segment_Mesh(); // edges only makes sense when referring to a segment mesh
    bool incident_segments_defined=segment_mesh->incident_elements!=0;if(!incident_segments_defined) segment_mesh->Initialize_Incident_Elements();
    for(int t=1;t<=elements.m;t++){
        int i,j,k;elements(t).Get(i,j,k);
        (*element_edges)(t).Set(segment_mesh->Segment(i,j),segment_mesh->Segment(j,k),segment_mesh->Segment(k,i));}
    if(!incident_segments_defined){delete segment_mesh->incident_elements;segment_mesh->incident_elements=0;}
}
//#####################################################################
// Function Initialize_Triangles_On_Edges
//#####################################################################
void TRIANGLE_MESH::
Initialize_Edge_Triangles()
{
    delete edge_triangles;
    if(!segment_mesh) Initialize_Segment_Mesh(); // edges only makes sense when referring to a segment mesh
    edge_triangles=new ARRAY<ARRAY<int> >(segment_mesh->elements.m);
    bool triangle_edges_defined=(element_edges!=0);if(!triangle_edges_defined) Initialize_Element_Edges();
    for(int t=1;t<=elements.m;t++){
        int e1,e2,e3;(*element_edges)(t).Get(e1,e2,e3);
        (*edge_triangles)(e1).Append(t);(*edge_triangles)(e2).Append(t);(*edge_triangles)(e3).Append(t);}
    if(!triangle_edges_defined){delete element_edges;element_edges=0;}
    for(int s=1;s<=segment_mesh->elements.m;s++){
        (*edge_triangles)(s).Compact();
        if((*edge_triangles)(s).m==2){
            int s1,s2,t1,t2,t3;segment_mesh->elements(s).Get(s1,s2);elements((*edge_triangles)(s)(1)).Get(t1,t2,t3);
            if(s2!=(s1==t1?t2:s1==t2?t3:t1))exchange((*edge_triangles)(s)(1),(*edge_triangles)(s)(2));}}
}
//#####################################################################
// Function Initialize_Boundary_Mesh
//#####################################################################
// The direction of each boundary segment matches that of the triangle it lives on
// TODO: Move into SIMPLEX_MESH somehow?
void TRIANGLE_MESH::
Initialize_Boundary_Mesh()
{
    delete boundary_mesh;
    boundary_mesh=new T_BOUNDARY_MESH();
    ARRAY< VECTOR<int,dimension> >& boundary_elements=boundary_mesh->elements;
    bool incident_elements_defined=incident_elements!=0;
    if(!incident_elements_defined) Initialize_Incident_Elements();
    for(int e=1;e<=elements.m;++e){
        const VECTOR<int,dimension+1>& element=elements(e);
        for(int i=1;i<=dimension+1;++i){
            VECTOR<int,dimension> face=element.Remove_Index(i);
            if(i==2) exchange(face.x,face.y); // ensure cyclic order
            const ARRAY<int>& incident_elements_to_face1=(*incident_elements)(face[1]);
            bool another_element_on_face=false;
            for(int j=1;j<=incident_elements_to_face1.m&&!another_element_on_face;++j)
                another_element_on_face=(incident_elements_to_face1(j)!=e&&Nodes_In_Simplex(face,incident_elements_to_face1(j)));
            if(!another_element_on_face) boundary_elements.Append(face);}}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
    boundary_elements.Compact();
    boundary_mesh->number_nodes=number_nodes;
}
//#####################################################################
// Function Initialize_Node_On_Boundary
//#####################################################################
// TODO: Move into SIMPLEX_MESH somehow?
void TRIANGLE_MESH::
Initialize_Node_On_Boundary()
{
    delete node_on_boundary;
    node_on_boundary=new ARRAY<bool>(number_nodes); // init'ed to false
    bool incident_elements_defined=incident_elements!=0;
    if(!incident_elements_defined) Initialize_Incident_Elements();
    for(int e=1;e<=elements.m;++e){
        const VECTOR<int,dimension+1>& element=elements(e);
        for(int i=1;i<=dimension+1;++i){
            VECTOR<int,dimension> face=element.Remove_Index(i);
            const ARRAY<int>& incident_elements_to_face1=(*incident_elements)(face[1]);
            bool another_element_on_face=false;
            for(int j=1;j<=incident_elements_to_face1.m&&!another_element_on_face;++j)
                another_element_on_face=(incident_elements_to_face1(j)!=e&&Nodes_In_Simplex(face,incident_elements_to_face1(j)));
            if(!another_element_on_face){INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,dimension>&> subset=node_on_boundary->Subset(face);ARRAYS_COMPUTATIONS::Fill(subset,true);}}}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Initialize_Boundary_Nodes
//#####################################################################
void TRIANGLE_MESH::
Initialize_Boundary_Nodes()
{
    delete boundary_nodes;boundary_nodes=new ARRAY<int>;
    bool boundary_mesh_defined=boundary_mesh!=0;if(!boundary_mesh_defined) Initialize_Boundary_Mesh();
    for(int t=1;t<=boundary_mesh->elements.m;t++){int i,j;boundary_mesh->elements(t).Get(i,j);boundary_nodes->Append(i);boundary_nodes->Append(j);}
    if(!boundary_mesh_defined){delete boundary_mesh;boundary_mesh=0;}

    Sort(*boundary_nodes);
    boundary_nodes->Prune_Duplicates();
}
//#####################################################################
// Function Non_Manifold_Nodes
//#####################################################################
void TRIANGLE_MESH::
Non_Manifold_Nodes(ARRAY<int>& node_list)
{
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    bool neighbor_nodes_defined=neighbor_nodes!=0;if(!neighbor_nodes) Initialize_Neighbor_Nodes();
    node_list.Remove_All();

    for(int i=1;i<=number_nodes;i++){
        if((*neighbor_nodes)(i).m != (*incident_elements)(i).m) node_list.Append(i);
        else if((*neighbor_nodes)(i).m > 0){
            ARRAY<int> ordered_neighbors;ordered_neighbors.Preallocate((*neighbor_nodes)(i).m);
            ordered_neighbors.Append((*neighbor_nodes)(i)(1));
            bool found_neighbor=true;
            while(found_neighbor){found_neighbor=false;
                int k=ordered_neighbors(ordered_neighbors.m); // looking for a new neighbor of node k
                for(int j=1;j<=(*incident_elements)(i).m;j++){
                    int a,b,c;elements((*incident_elements)(i)(j)).Get(a,b,c);
                    if(a == k || b == k || c == k){int index;
                        if(a != i && !ordered_neighbors.Find(a,index)){ordered_neighbors.Append(a);found_neighbor=true;break;}
                        else if(b != i && !ordered_neighbors.Find(b,index)){ordered_neighbors.Append(b);found_neighbor=true;break;}
                        else if(c != i && !ordered_neighbors.Find(c,index)){ordered_neighbors.Append(c);found_neighbor=true;break;}}}}
            if(ordered_neighbors.m < (*neighbor_nodes)(i).m) node_list.Append(i);}}

    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
    if(!neighbor_nodes_defined){delete neighbor_nodes;neighbor_nodes=0;}
}
//#####################################################################
// Function Mark_Face_Connected_Component_Incident_On_A_Node
//#####################################################################
// implements a depth first search to find marked connected component
void TRIANGLE_MESH::
Mark_Edge_Connected_Component_Incident_On_A_Node(const int node,const int triangle_index_in_incident_elements,ARRAY<bool>& marked) const
{
    assert(incident_elements); // without this triangle_index_in_incident_elements makes no sense
    assert(adjacent_elements); // too expensive to do this every time, since this is a recursive function
    int triangle=(*incident_elements)(node)(triangle_index_in_incident_elements);
    if(marked.m != (*incident_elements)(node).m) marked.Resize((*incident_elements)(node).m);
    marked(triangle_index_in_incident_elements)=true;
    for(int i=1;i<=(*adjacent_elements)(triangle).m;i++){
        int adjacent_triangle=(*adjacent_elements)(triangle)(i);
        if(Node_In_Triangle(node,adjacent_triangle)){int index;
            if((*incident_elements)(node).Find(adjacent_triangle,index) && !marked(index)) Mark_Edge_Connected_Component_Incident_On_A_Node(node,index,marked);}}
}
//#####################################################################
// Function Adjacent_Triangle
//#####################################################################
int TRIANGLE_MESH::
Adjacent_Triangle(const int triangle,const int node1,const int node2) const
{
    assert(Node_In_Triangle(node1,triangle) && Node_In_Triangle(node2,triangle) && adjacent_elements);
    for(int a=1;a<=(*adjacent_elements)(triangle).m;a++){
        int t=(*adjacent_elements)(triangle)(a);if(Node_In_Triangle(node1,t) && Node_In_Triangle(node2,t))return t;}
    return 0;
}
//#####################################################################
// Function Triangles_On_Edge
//#####################################################################
int TRIANGLE_MESH::
Triangles_On_Edge(const int node1,const int node2,ARRAY<int>* triangles_on_edge) const
{
    assert(incident_elements);
    int triangle_count=0;
    for(int t=1;t<=(*incident_elements)(node1).m;t++){
        int triangle=(*incident_elements)(node1)(t);
        if(Node_In_Triangle(node2,triangle)){triangle_count++;if(triangles_on_edge) (*triangles_on_edge).Append(triangle);}}
    return triangle_count;
}
//#####################################################################
// Function Triangle_On_Edge
//#####################################################################
bool TRIANGLE_MESH::
Triangle_On_Edge(const int node1,const int node2) const
{
    assert(incident_elements);
    for(int t=1;t<=(*incident_elements)(node1).m;t++) if(Node_In_Triangle(node2,(*incident_elements)(node1)(t)))return true;
    return false;
}
//#####################################################################
// Function Triangles_Across_Edge
//#####################################################################
int TRIANGLE_MESH::
Triangles_Across_Edge(const int triangle,const int node1,const int node2,ARRAY<int>& triangles_across_edge) const
{
    assert(incident_elements);
    for(int k=1;k<=(*incident_elements)(node1).m;k++){
        int triangle2=(*incident_elements)(node1)(k); // triangle in question
        if(triangle2 != triangle && Node_In_Triangle(node2,triangle2))triangles_across_edge.Append(triangle2);}
    return triangles_across_edge.m;
}
//#####################################################################
// Function Make_Orientations_Consistent
//#####################################################################
void TRIANGLE_MESH::
Make_Orientations_Consistent()
{
    bool adjacent_elements_defined=adjacent_elements!=0;if(!adjacent_elements_defined)Initialize_Adjacent_Elements();
    ARRAY<int> orientation(elements.m),component(1000);
    for(int c=1;c<=elements.m;c++)if(!orientation(c)){
        int positive=0,negative=0;
        component.Remove_All();component.Append(c);orientation(c)=1;
        for(int a=1;a<=component.m;a++){
            int t=component(a),i,j,k;elements(t).Get(i,j,k);
            if(orientation(c)>0)positive++;else negative++;
            for(int b=1;b<=(*adjacent_elements)(t).m;b++){
                int t2=(*adjacent_elements)(t)(b);if(orientation(t2))continue;
                int i2,j2,k2;elements(t2).Get(i2,j2,k2);if(i!=i2&&i!=j2&&i!=k2)cyclic_shift(i,j,k);
                if(i!=i2)cyclic_shift(i2,j2,k2);if(i!=i2)cyclic_shift(i2,j2,k2);
                component.Append(t2);orientation(t2)=(j==k2||k==j2)?orientation(t):-orientation(t);}}
        int correct=positive>negative?1:-1;
        for(int a=1;a<=component.m;a++)if(orientation(component(a))!=correct)
            exchange(elements(component(a))(2),elements(component(a))(3));}
    if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}
}
//#####################################################################
// Function Make_Orientations_Consistent
//#####################################################################
bool TRIANGLE_MESH::
Orientations_Consistent()
{
    bool adjacent_elements_defined=adjacent_elements!=0;if(!adjacent_elements_defined)Initialize_Adjacent_Elements();
    for(int t=1;t<=elements.m;t++){
        int i,j,k;elements(t).Get(i,j,k);
        for(int a=1;a<=(*adjacent_elements)(t).m;a++){
            int t2=(*adjacent_elements)(t)(a);if(t2<t)continue;
            int i2,j2,k2;elements(t2).Get(i2,j2,k2);if(i!=i2&&i!=j2&&i!=k2)cyclic_shift(i,j,k);
            if(i!=i2)cyclic_shift(i2,j2,k2);if(i!=i2)cyclic_shift(i2,j2,k2);
            if(j==j2||k==k2){if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}return false;}}}
    if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}return true;
}
//#####################################################################
// Function Identify_Connected_Components
//#####################################################################
// Labels each connected components with a unique id
void TRIANGLE_MESH::
Identify_Connected_Components(ARRAY<int>& label)
{
    bool adjacent_elements_defined=adjacent_elements!=0;if(!adjacent_elements_defined)Initialize_Adjacent_Elements();
    label.Resize(elements.m);ARRAYS_COMPUTATIONS::Fill(label,0);
    int id=0;for(int i=1;i<=label.m;i++) if(!label(i)) Label_Connected_Component_With_ID(label,i,++id);
    if(!adjacent_elements_defined){delete adjacent_elements;adjacent_elements=0;}
}
//#####################################################################
// Function Label_Connected_Component_With_ID
//#####################################################################
// Labels connected components to a given triangle with a unique id
void TRIANGLE_MESH::
Label_Connected_Component_With_ID(ARRAY<int>& label,const int triangle,const int id) const
{
    assert(adjacent_elements);assert(label.m==elements.m);
    label(triangle)=id;
    for(int i=1;i<=(*adjacent_elements)(triangle).m;i++){
        int adjacent_triangle=(*adjacent_elements)(triangle)(i);
        if(!label(adjacent_triangle)) Label_Connected_Component_With_ID(label,adjacent_triangle,id);}
}
//#####################################################################
// Function Assert_Consistent
//#####################################################################
bool TRIANGLE_MESH::
Assert_Consistent() const
{
    if(topologically_sorted_neighbor_nodes) assert(topologically_sorted_neighbor_nodes->m==number_nodes);
    if(topologically_sorted_incident_elements) assert(topologically_sorted_incident_elements->m==number_nodes);
    if(segment_mesh){assert(segment_mesh->number_nodes==number_nodes);assert(segment_mesh->Assert_Consistent());}
    if(element_edges){assert(segment_mesh);assert(element_edges->m==elements.m);}
    if(edge_triangles){assert(segment_mesh);assert(edge_triangles->m==segment_mesh->elements.m);}
    if(boundary_mesh){assert(boundary_mesh->number_nodes==number_nodes);assert(boundary_mesh->Assert_Consistent());}
    if(node_on_boundary) assert(node_on_boundary->m==number_nodes);
    return SIMPLEX_MESH<2>::Assert_Consistent();
}
//#####################################################################
// Function Set_Number_Nodes
//#####################################################################
void TRIANGLE_MESH::
Set_Number_Nodes(const int number_nodes_input)
{
    SIMPLEX_MESH<2>::Set_Number_Nodes(number_nodes_input);
    if(segment_mesh) segment_mesh->Set_Number_Nodes(number_nodes);
    if(boundary_mesh) boundary_mesh->Set_Number_Nodes(number_nodes);
    if(node_on_boundary) node_on_boundary->Resize(number_nodes);
}
//#####################################################################
