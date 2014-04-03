//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLEX_MESH
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#endif
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
#include <climits>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<int d> SIMPLEX_MESH<d>::
SIMPLEX_MESH()
    :number_nodes(0),neighbor_nodes(0),incident_elements(0),adjacent_elements(0),neighbor_elements(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<int d> SIMPLEX_MESH<d>::
SIMPLEX_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,d+1> >& simplex_list)
    :neighbor_nodes(0),incident_elements(0),adjacent_elements(0),neighbor_elements(0)
{
    Initialize_Mesh(number_nodes_input,simplex_list);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> SIMPLEX_MESH<d>::
SIMPLEX_MESH(const SIMPLEX_MESH& mesh)
    :neighbor_nodes(0),incident_elements(0),adjacent_elements(0),neighbor_elements(0)
{
    Initialize_Mesh(mesh);
}
//#####################################################################
// Destructor
//#####################################################################
template<int d> SIMPLEX_MESH<d>::
~SIMPLEX_MESH()
{
    delete neighbor_nodes;delete incident_elements;delete adjacent_elements;delete neighbor_elements;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Clean_Memory()
{
    elements.Clean_Memory();Delete_Auxiliary_Structures();
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Delete_Auxiliary_Structures()
{
    delete neighbor_nodes;neighbor_nodes=0;delete incident_elements;incident_elements=0;delete adjacent_elements;adjacent_elements=0;delete neighbor_elements;neighbor_elements=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Refresh_Auxiliary_Structures()
{
    if(neighbor_nodes) Initialize_Neighbor_Nodes();if(incident_elements) Initialize_Incident_Elements();
    if(adjacent_elements) Initialize_Adjacent_Elements();if(neighbor_elements) Initialize_Neighbor_Elements();
}
//#####################################################################
// Function Simplex
//#####################################################################
// requires incident_elements
// returns the simplex containing these nodes - returns 0 if no simplex contains these nodes
template<int d> int SIMPLEX_MESH<d>::
Simplex(const VECTOR<int,d+1>& nodes) const
{
    assert(incident_elements);
    // find shortest list of elements
    int short_list=nodes[1];VECTOR<int,d> check=nodes.Remove_Index(1);
    for(int i=1;i<=d;i++)if((*incident_elements)(check[i]).m < (*incident_elements)(short_list).m) exchange(short_list,check[i]);
    // search short list for other nodes
    for(int k=1;k<=(*incident_elements)(short_list).m;k++){int t=(*incident_elements)(short_list)(k);
        if(Nodes_In_Simplex(check,t)) return t;}
    return 0;
}
//#####################################################################
// Function Initialize_Neighbor_Nodes
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Initialize_Neighbor_Nodes()
{
    delete neighbor_nodes;neighbor_nodes=new ARRAY<ARRAY<int> >(number_nodes);
    for(int t=1;t<=elements.m;t++){VECTOR<int,d+1>& element=elements(t);
        for(int i=1;i<=element.m;i++)for(int j=1;j<=element.m;j++)if(i!=j)
            (*neighbor_nodes)(element[i]).Append_Unique(element[j]);}
    for(int p=1;p<=number_nodes;p++) (*neighbor_nodes)(p).Compact(); // remove extra space
}
//#####################################################################
// Function Initialize_Incident_Elements
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Initialize_Incident_Elements()
{
    delete incident_elements;incident_elements=new ARRAY<ARRAY<int> >(number_nodes);
    for(int t=1;t<=elements.m;t++){VECTOR<int,d+1>& element=elements(t); // for each element, put it on each of its nodes' lists of simplices
        for(int i=1;i<=element.m;i++) (*incident_elements)(element[i]).Append(t);}
    for(int p=1;p<=number_nodes;p++) (*incident_elements)(p).Compact(); // remove extra space
}
//#####################################################################
// Function Initialize_Adjacent_Elements
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Initialize_Adjacent_Elements()
{
    delete adjacent_elements;adjacent_elements=new ARRAY<ARRAY<int> >(elements.m);
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int t=1;t<=elements.m;t++){VECTOR<int,d+1>& element=elements(t); // for each simplex, make the list of adjacent simplices
        for(int i=1;i<=element.m;i++) Find_And_Append_Adjacent_Elements(t,element.Remove_Index(i));}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
    for(int t=1;t<=elements.m;t++) (*adjacent_elements)(t).Compact(); // remove extra space
}
//#####################################################################
// Function Find_And_Append_Adjacent_Elements
//#####################################################################
// find adjacent simplices that contains face and append them to the adjacency list
template<int d> void SIMPLEX_MESH<d>::
Find_And_Append_Adjacent_Elements(const int simplex,const VECTOR<int,d>& face)
{
    int first_node=face[1];VECTOR<int,d-1> other_nodes=face.Remove_Index(1);
    for(int t=1;t<=(*incident_elements)(first_node).m;t++){
        int simplex2=(*incident_elements)(first_node)(t);
        if(simplex2!=simplex && Nodes_In_Simplex(other_nodes,simplex2))
            (*adjacent_elements)(simplex).Append_Unique(simplex2);}
}
template<> void SIMPLEX_MESH<0>:: // specialize segment case since other_nodes would be 0 size
Find_And_Append_Adjacent_Elements(const int simplex,const VECTOR<int,0>& face)
{return;}
//#####################################################################
// Function Initialize_Neighbor_Elements
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Initialize_Neighbor_Elements()
{
    STATIC_ASSERT(d<4);static const int expected_neighbors=d==1?4:d==2?20:d==3?50:0;
    delete neighbor_elements;neighbor_elements=new ARRAY<ARRAY<int> >(elements.m);
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int t=1;t<=elements.m;t++){VECTOR<int,d+1>& element=elements(t); // for each simplex, make the list of neighbor simplices
        ARRAY<int>& neighbors=(*neighbor_elements)(t);neighbors.Preallocate(expected_neighbors);
        neighbors.Append_Elements((*incident_elements)(element[1]));
        for(int i=2;i<=element.m;i++) neighbors.Append_Unique_Elements((*incident_elements)(element[i]));
        neighbors.Remove_Index_Lazy(neighbors.Find(t));neighbors.Compact();}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Add_Element_If_Not_Already_There
//#####################################################################
// will ONLY incrementally update incident_elements, assumes number_nodes is already set to accommodate the new element
template<int d> int SIMPLEX_MESH<d>::
Add_Element_If_Not_Already_There(const VECTOR<int,d+1>& nodes)
{
    if(int element=Simplex(nodes)) return element;
    elements.Append(nodes);
    for(int i=1;i<=d+1;i++) (*incident_elements)(nodes[i]).Append(elements.m);
    return elements.m;
}
//#####################################################################
// Function Delete_Elements_With_Missing_Nodes
//#####################################################################
template<int d> int SIMPLEX_MESH<d>::
Delete_Elements_With_Missing_Nodes()
{
    int m_save=elements.m;
    for(int t=elements.m;t>=1;t--)if(elements(t).Contains(0)) elements.Remove_Index_Lazy(t);
    elements.Compact();Refresh_Auxiliary_Structures();
    return m_save-elements.m;
}
//#####################################################################
// Function Delete_Sorted_Elements
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Delete_Sorted_Elements(const ARRAY<int>& deletion_list)
{
    elements.Remove_Sorted_Indices_Lazy(deletion_list);
    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Delete_Sorted_Elements
// index_map. Maps from old element index to new element index. (0 if deleted, no entry if unchanged)
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Delete_Sorted_Elements(const ARRAY<int>& deletion_list,HASHTABLE<int,int>& index_map) 
{
    int curr=0;
    HASHTABLE<int,int> reverse_index_map;
    ARRAY<int> unique_deletion_list;
    for(int k=deletion_list.m;k>=1;k--) if(deletion_list(k)!=curr){
        curr=deletion_list(k);int previous_index=elements.m;
        unique_deletion_list.Append(curr);
        if(curr!=elements.m){
            if(reverse_index_map.Get(elements.m,previous_index)) reverse_index_map.Delete(elements.m);
            reverse_index_map.Insert(curr,previous_index);}
        elements.Remove_Index_Lazy(curr);}
    for(HASHTABLE_ITERATOR<int,int> iter(reverse_index_map);iter.Valid();iter.Next()) index_map.Insert(iter.Data(),iter.Key());
    for(int i=1;i<=unique_deletion_list.m;i++) index_map.Insert(unique_deletion_list(i),0);
    elements.Compact();Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Delete_Elements
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Delete_Elements(ARRAY<int> deletion_list)
{
    Sort(deletion_list);Delete_Sorted_Elements(deletion_list);
}
//#####################################################################
// Function Number_Of_Nodes_With_Minimum_Valence
//#####################################################################
template<int d> int SIMPLEX_MESH<d>::
Number_Of_Nodes_With_Minimum_Valence()
{
    bool neighbor_nodes_defined=neighbor_nodes!=0;if(!neighbor_nodes_defined) Initialize_Neighbor_Nodes();
    int number_of_nodes_with_minimum_valence=0,minimum=Minimum_Valence();
    for(int t=1;t<=number_nodes;t++) if(minimum == (*neighbor_nodes)(t).m) number_of_nodes_with_minimum_valence++;
    if(!neighbor_nodes_defined){delete neighbor_nodes;neighbor_nodes=0;}
    return number_of_nodes_with_minimum_valence;
}
//#####################################################################
// Function Minimum_Valence
//#####################################################################
template<int d> int SIMPLEX_MESH<d>::
Minimum_Valence(int* index)
{
    bool neighbor_nodes_defined=neighbor_nodes!=0;if(!neighbor_nodes_defined) Initialize_Neighbor_Nodes();
    int minimum=INT_MAX,minimum_index=0;
    for(int p=1;p<=number_nodes;p++)if((*neighbor_nodes)(p).m && (*neighbor_nodes)(p).m < minimum){
        minimum=(*neighbor_nodes)(p).m;minimum_index=p;}
    if(index) *index=minimum_index;
    if(!neighbor_nodes_defined){delete neighbor_nodes;neighbor_nodes=0;}
    return minimum;
}
//#####################################################################
// Function Maximum_Valence
//#####################################################################
template<int d> int SIMPLEX_MESH<d>::
Maximum_Valence(int* index)
{
    bool neighbor_nodes_defined=neighbor_nodes!=0;if(!neighbor_nodes_defined) Initialize_Neighbor_Nodes();
    int maximum=0,maximum_index=0;
    for(int p=1;p<=number_nodes;p++)if((*neighbor_nodes)(p).m > maximum){
        maximum=(*neighbor_nodes)(p).m;maximum_index=p;}
    if(index) *index=maximum_index;
    if(!neighbor_nodes_defined){delete neighbor_nodes;neighbor_nodes=0;}
    return maximum;
}
//#####################################################################
// Function Update_Adjacent_Elements_From_Incident_Elements
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Update_Adjacent_Elements_From_Incident_Elements(const int node)
{
    if(!adjacent_elements){Initialize_Adjacent_Elements();return;}
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int a=1;a<=(*incident_elements)(node).m;a++){
        int t=(*incident_elements)(node)(a);VECTOR<int,d+1>& element=elements(t);
        (*adjacent_elements)(t).Remove_All();
        for(int i=1;i<=element.m;i++) Find_And_Append_Adjacent_Elements(t,element.Remove_Index(i));
        (*adjacent_elements)(t).Compact();}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Update_Neighbor_Nodes_From_Incident_Elements
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Update_Neighbor_Nodes_From_Incident_Elements(const int node)
{
    if(!neighbor_nodes){Initialize_Neighbor_Nodes();return;}
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    (*neighbor_nodes)(node).Remove_All();
    for(int t=1;t<=(*incident_elements)(node).m;t++){
        VECTOR<int,d+1>& element=elements((*incident_elements)(node)(t));assert(element.Contains(node));
        for(int i=1;i<=element.m;i++)if(element[i]!=node) (*neighbor_nodes)(node).Append_Unique(element[i]);}
    (*neighbor_nodes)(node).Compact();
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Assert_Consistent
//#####################################################################
template<int d> bool SIMPLEX_MESH<d>::
Assert_Consistent() const
{
    if(neighbor_nodes) assert(neighbor_nodes->m==number_nodes);
    if(incident_elements) assert(incident_elements->m==number_nodes);
    if(adjacent_elements) assert(adjacent_elements->m==elements.m);
    if(neighbor_elements) assert(neighbor_elements->m==elements.m);
    return true;
}
//#####################################################################
// Function Set_Number_Nodes
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Set_Number_Nodes(const int number_nodes_input)
{
    assert(Assert_Consistent());
    number_nodes=number_nodes_input;
    if(neighbor_nodes) neighbor_nodes->Resize(number_nodes);
    if(incident_elements) incident_elements->Resize(number_nodes);
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<int d> void SIMPLEX_MESH<d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int t=1;t<=elements.m;t++) for(int i=1;i<=d;i++) for(int j=i+1;j<=d+1;j++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(elements(t)[i],elements(t)[j]));
}
//#####################################################################
// Function Mark_Nodes_Referenced
//#####################################################################
template<int d> template<class T> void SIMPLEX_MESH<d>::
Mark_Nodes_Referenced(ARRAY<T>& marks,const T& mark) const
{
    for(int e=1;e<=elements.m;e++){INDIRECT_ARRAY<ARRAY<T>,VECTOR<int,d+1>&> subset=marks.Subset(elements(e));ARRAYS_COMPUTATIONS::Fill(subset,mark);}
}
//#####################################################################
// Function Simplices_On_Subsimplex
//#####################################################################
template<int d> template<int d2> void SIMPLEX_MESH<d>::
Simplices_On_Subsimplex(const VECTOR<int,d2>& subsimplex_nodes,ARRAY<int>& simplices_on_subsimplex) const
{
    assert(incident_elements);
    const ARRAY<int>& incident=(*incident_elements)(subsimplex_nodes[1]);
    VECTOR<int,d2-1> other_nodes=subsimplex_nodes.Remove_Index(1);
    for(int i=1;i<=incident.m;i++){
        int simplex=incident(i);
        if(Nodes_In_Simplex(other_nodes,simplex)) simplices_on_subsimplex.Append(simplex);}
}
//#####################################################################
template class SIMPLEX_MESH<0>;
template class SIMPLEX_MESH<1>;
template class SIMPLEX_MESH<2>;
template class SIMPLEX_MESH<3>;
template void SIMPLEX_MESH<2>::Simplices_On_Subsimplex(const VECTOR<int,2>&,ARRAY<int>&) const;
template void SIMPLEX_MESH<3>::Simplices_On_Subsimplex(const VECTOR<int,2>&,ARRAY<int>&) const;
template void SIMPLEX_MESH<3>::Simplices_On_Subsimplex(const VECTOR<int,3>&,ARRAY<int>&) const;
template void SIMPLEX_MESH<0>::Mark_Nodes_Referenced(ARRAY<int>&,const int&) const;
template void SIMPLEX_MESH<1>::Mark_Nodes_Referenced(ARRAY<int>&,const int&) const;
template void SIMPLEX_MESH<2>::Mark_Nodes_Referenced(ARRAY<int>&,const int&) const;
template void SIMPLEX_MESH<3>::Mark_Nodes_Referenced(ARRAY<int>&,const int&) const;
template void SIMPLEX_MESH<3>::Mark_Nodes_Referenced(ARRAY<bool>&,const bool&) const;
}
