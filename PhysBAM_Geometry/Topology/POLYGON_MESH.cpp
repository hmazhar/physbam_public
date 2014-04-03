//#####################################################################
// Copyright 2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGON_MESH
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Topology/POLYGON_MESH.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
POLYGON_MESH::
POLYGON_MESH()
    :number_nodes(0),segment_mesh(0),element_oriented_edges(0),edge_elements(0)
{}
//#####################################################################
// Constructor
//#####################################################################
POLYGON_MESH::
POLYGON_MESH(const int number_nodes_input,const ARRAY<ARRAY<ARRAY<int> > >& polygon_list)
    :segment_mesh(0),element_oriented_edges(0),edge_elements(0)
{
    Initialize_Mesh(number_nodes_input,polygon_list);
}
//#####################################################################
// Constructor
//#####################################################################
POLYGON_MESH::
POLYGON_MESH(const POLYGON_MESH& mesh)
    :segment_mesh(0),element_oriented_edges(0),edge_elements(0)
{
    Initialize_Mesh(mesh);
}
//#####################################################################
// Destructor
//#####################################################################
POLYGON_MESH::
~POLYGON_MESH()
{
    delete segment_mesh;delete element_oriented_edges;delete edge_elements;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
void POLYGON_MESH::
Clean_Memory()
{
    elements.Clean_Memory();Delete_Auxiliary_Structures();
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
void POLYGON_MESH::
Delete_Auxiliary_Structures()
{
    delete segment_mesh;segment_mesh=0;delete element_oriented_edges;element_oriented_edges=0;delete edge_elements;edge_elements=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures
//#####################################################################
void POLYGON_MESH::
Refresh_Auxiliary_Structures()
{
    if(segment_mesh){
        bool incident_elements_is_defined=segment_mesh->incident_elements!=0;
        Initialize_Segment_Mesh();
        if(incident_elements_is_defined) segment_mesh->Initialize_Incident_Elements();}
    if(element_oriented_edges) Initialize_Element_Oriented_Edges();
    if(edge_elements) Initialize_Edge_Elements();
}
//#####################################################################
// Function Assert_Consistent
//#####################################################################
bool POLYGON_MESH::
Assert_Consistent() const
{
    if(segment_mesh) assert(segment_mesh->Assert_Consistent());
    return true;
}
//#####################################################################
// Function Set_Number_Nodes
//#####################################################################
void POLYGON_MESH::
Set_Number_Nodes(const int number_nodes_input)
{
    assert(Assert_Consistent());
    if(number_nodes_input<number_nodes) PHYSBAM_FATAL_ERROR();
    number_nodes=number_nodes_input;
    if(segment_mesh) segment_mesh->Set_Number_Nodes(number_nodes);
}
//#####################################################################
// Function Initialize_Segment_Mesh
//#####################################################################
void POLYGON_MESH::
Initialize_Segment_Mesh()
{
    delete segment_mesh;segment_mesh=new SEGMENT_MESH;segment_mesh->number_nodes=number_nodes;
    HASHTABLE<VECTOR<int,2> > segment_list;
    for(int e=1;e<=elements.m;e++) for(int c=1;c<=elements(e).m;c++){
        ARRAY<int>& component=elements(e)(c);
        for(int p=1;p<=component.m;p++){
            VECTOR<int,2> sorted_segment=VECTOR<int,2>(component(p),component(p%component.m+1)).Sorted();
            if(segment_list.Set(sorted_segment)) segment_mesh->elements.Append(sorted_segment);}}
}
//#####################################################################
// Function Initialize_Element_Oriented_Edges
//#####################################################################
void POLYGON_MESH::
Initialize_Element_Oriented_Edges()
{
    if(!segment_mesh) PHYSBAM_FATAL_ERROR();
    delete element_oriented_edges;element_oriented_edges=new ARRAY<ARRAY<ARRAY<PAIR<int,bool> > > >(elements.m);
    for(int e=1;e<=elements.m;e++){
        (*element_oriented_edges)(e).Resize(elements(e).m);
        for(int c=1;c<=elements(e).m;c++){
            ARRAY<int>& component=elements(e)(c);
            for(int p=1;p<=component.m;p++){
                VECTOR<int,2> oriented_segment(component(p),component(p%component.m+1));
                int s=segment_mesh->Simplex(oriented_segment);if(!s) PHYSBAM_FATAL_ERROR();
                if(oriented_segment==segment_mesh->elements(s)) (*element_oriented_edges)(e)(c).Append(PAIR<int,bool>(s,true));
                else (*element_oriented_edges)(e)(c).Append(PAIR<int,bool>(s,false));}}}
}
//#####################################################################
// Function Initialize_Edge_Elements
//#####################################################################
void POLYGON_MESH::
Initialize_Edge_Elements()
{
    if(!segment_mesh) PHYSBAM_FATAL_ERROR();
    delete edge_elements;edge_elements=new ARRAY<ARRAY<int> >(segment_mesh->elements.m);
    OPERATION_HASH<> hash(segment_mesh->elements.m);
    for(int e=1;e<=elements.m;e++) for(int c=1;c<=elements(e).m;c++){
        ARRAY<int>& component=elements(e)(c);
        for(int p=1;p<=component.m;p++){
            int s=segment_mesh->Segment(component(p),component(p%component.m+1));if(!s) PHYSBAM_FATAL_ERROR();
            if(!hash.Is_Marked_Current(s)){(*edge_elements)(s).Append(e);hash.Mark(s);}}
        hash.Next_Operation();}
}
//#####################################################################
// Function Elements_On_Edge
//#####################################################################
int POLYGON_MESH::
Elements_On_Edge(const int node1,const int node2,ARRAY<int>* elements_on_edge) const
{
    if(!segment_mesh || !edge_elements) PHYSBAM_FATAL_ERROR();
    int s=segment_mesh->Segment(node1,node2);
    if(!s){if(elements_on_edge) elements_on_edge->Remove_All();return 0;}
    if(elements_on_edge) *elements_on_edge=(*edge_elements)(s);
    return (*edge_elements)(s).m;
}
//#####################################################################
// Function Oriented_Edge_In_Element
//#####################################################################
bool POLYGON_MESH::
Oriented_Edge_In_Element(const int node1,const int node2,const int element) const
{
    if(!segment_mesh || !element_oriented_edges) PHYSBAM_FATAL_ERROR();
    int s=segment_mesh->Segment(node1,node2);if(!s) return false;
    PAIR<int,bool> oriented_edge(s,segment_mesh->elements(s)==VECTOR<int,2>(node1,node2));
    for(int c=1;c<=(*element_oriented_edges)(element).m;c++)
        if((*element_oriented_edges)(element)(c).Contains(oriented_edge)) return true;
    return false;
}
//#####################################################################
// Function Elements_On_Oriented_Edge
//#####################################################################
int POLYGON_MESH::
Elements_On_Oriented_Edge(const int node1,const int node2,ARRAY<int>* elements_on_oriented_edge) const
{
    if(!segment_mesh || !edge_elements || !element_oriented_edges) PHYSBAM_FATAL_ERROR();
    int s=segment_mesh->Segment(node1,node2);
    if(!s){if(elements_on_oriented_edge) elements_on_oriented_edge->Remove_All();return 0;}
    ARRAY<int>& candidate_elements_on_edge=(*edge_elements)(s);
    int elements_found=0;if(elements_on_oriented_edge) elements_on_oriented_edge->Remove_All();
    for(int i=1;i<=candidate_elements_on_edge.m;i++) if(Oriented_Edge_In_Element(node1,node2,candidate_elements_on_edge(i))){
        elements_found++;if(elements_on_oriented_edge) elements_on_oriented_edge->Append(candidate_elements_on_edge(i));}
    return elements_found;
}
//#####################################################################
// Function Opposite_Oriented_Element
//#####################################################################
// Assumes that polygons are connected (but not necessarily simply connected)
int POLYGON_MESH::
Opposite_Oriented_Element(const int element) const
{
    if(!segment_mesh || !element_oriented_edges) PHYSBAM_FATAL_ERROR();
    ARRAY<int> candidate_elements=(*edge_elements)((*element_oriented_edges)(element)(1)(1).x);
    for(int c=1;c<=(*element_oriented_edges)(element).m;c++){
        for(int s=(c==1)?2:1;s<=(*element_oriented_edges)(element)(c).m;s++){
            PAIR<int,bool>& oriented_edge=(*element_oriented_edges)(element)(c)(s);
            int node1,node2;segment_mesh->elements(oriented_edge.x).Get(node1,node2);
            if(oriented_edge.y) exchange(node1,node2);
            for(int i=candidate_elements.m;i>=1;i--) if(!Oriented_Edge_In_Element(node1,node2,candidate_elements(i))) candidate_elements.Remove_Index_Lazy(i);}}
    for(int candidate=candidate_elements.m;candidate>=1;candidate--){
        for(int c=1;c<=(*element_oriented_edges)(candidate_elements(candidate)).m;c++)
            for(int s=1;s<=(*element_oriented_edges)(candidate_elements(candidate))(c).m;s++){
                PAIR<int,bool>& oriented_edge=(*element_oriented_edges)(candidate_elements(candidate))(c)(s);
                int node1,node2;segment_mesh->elements(oriented_edge.x).Get(node1,node2);
                if(oriented_edge.y) exchange(node1,node2);
                if(!Oriented_Edge_In_Element(node1,node2,element)){candidate_elements.Remove_Index_Lazy(candidate);goto Next_Candidate;}}
        Next_Candidate:;}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(candidate_elements.m>1){ 
        LOG::cout << "Input: " << elements(element) << std::endl;
        for(int i=1;i<=candidate_elements.m;i++)
            LOG::cout << "Candidate: " << elements(candidate_elements(i)) << std::endl;}
#endif
    if(candidate_elements.m>1) PHYSBAM_FATAL_ERROR();
    if(candidate_elements.m) return candidate_elements(1); else return 0;
}
//#####################################################################
// Function Split_Polygon_Edge
//#####################################################################
void POLYGON_MESH::
Split_Polygon_Edge(const int node1,const int node2,const int new_node)
{
    assert(segment_mesh->Simplex(VECTOR<int,2>(node1,node2)));
    // Incremental Update
    if(!segment_mesh || !segment_mesh->incident_elements || !element_oriented_edges || !edge_elements) PHYSBAM_FATAL_ERROR();
    if(number_nodes<new_node) Set_Number_Nodes(new_node);
    // compute elements on edge to be used later
    ARRAY<int> elements_on_unoriented_edge;Elements_On_Edge(node1,node2,&elements_on_unoriented_edge);
    // first edge replaces parent edge, second edge is added at the end
    int first_edge_index=segment_mesh->Segment(node1,node2),second_edge_index=segment_mesh->elements.Append(VECTOR<int,2>());
    VECTOR<int,2> &first_edge=segment_mesh->elements(first_edge_index),&second_edge=segment_mesh->elements(second_edge_index);
    if(first_edge==VECTOR<int,2>(node1,node2)){first_edge=VECTOR<int,2>(node1,new_node);second_edge=VECTOR<int,2>(new_node,node2);}
    else{first_edge=VECTOR<int,2>(node2,new_node);second_edge=VECTOR<int,2>(new_node,node1);}
    int index=(*segment_mesh->incident_elements)(second_edge.y).Find(first_edge_index);if(!index) PHYSBAM_FATAL_ERROR();
    (*segment_mesh->incident_elements)(second_edge.y)(index)=second_edge_index;
    (*segment_mesh->incident_elements)(new_node).Append(first_edge_index);(*segment_mesh->incident_elements)(new_node).Append(second_edge_index);
    // update elements and edge_elements
    ARRAY<int> elements_on_edge=(*edge_elements)(first_edge_index);edge_elements->Append(elements_on_edge);
    // update elements_on_unoriented_edge
    for(int e=1;e<=elements_on_unoriented_edge.m;e++){int element=elements_on_unoriented_edge(e);
        for(int c=1;c<=(*element_oriented_edges)(element).m;c++){
            for(int i=1;i<=(*element_oriented_edges)(element)(c).m;i++){
                PAIR<int,bool>& oriented_edge=(*element_oriented_edges)(element)(c)(i);
                if(oriented_edge.x==first_edge_index){
                    if(oriented_edge.y){
                        (*element_oriented_edges)(element)(c).Insert(PAIR<int,bool>(second_edge_index,true),++i);}
                    else{
                        (*element_oriented_edges)(element)(c).Insert(PAIR<int,bool>(second_edge_index,false),i++);}
                    elements(element)(c).Insert(new_node,i);}}}}
}
//#####################################################################
}
