//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRON_MESH
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Math_Tools/exchange_sort.h>
#include <PhysBAM_Geometry/Topology/HEXAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
using namespace PhysBAM;
//#####################################################################
const int HEXAHEDRON_MESH::face_indices[6][4]={{1,2,4,3},{5,7,8,6},{3,4,8,7},{1,5,6,2},{2,6,8,4},{1,3,7,5}};
const int HEXAHEDRON_MESH::edge_indices[12][2]={{1,2},{2,4},{4,3},{3,1},{1,5},{2,6},{4,8},{3,7},{5,6},{6,8},{8,7},{7,5}};
//#####################################################################
// Constructor
//#####################################################################
HEXAHEDRON_MESH::
HEXAHEDRON_MESH()
    :number_nodes(0),incident_elements(0),adjacent_elements(0),faces(0),node_on_boundary(0),boundary_nodes(0),face_hexahedrons(0),boundary_mesh(0)
{}
//#####################################################################
// Constructor
//#####################################################################
HEXAHEDRON_MESH::
HEXAHEDRON_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,8> >& hexahedron_list)
    :incident_elements(0),adjacent_elements(0),faces(0),node_on_boundary(0),face_hexahedrons(0),boundary_mesh(0)
{
    Initialize_Mesh(number_nodes_input,hexahedron_list);
}
//#####################################################################
// Destructor
//#####################################################################
HEXAHEDRON_MESH::
~HEXAHEDRON_MESH()
{
    delete incident_elements;delete adjacent_elements;delete faces;delete node_on_boundary;delete boundary_mesh;delete face_hexahedrons;
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
void HEXAHEDRON_MESH::
Delete_Auxiliary_Structures()
{
    delete incident_elements;incident_elements=0;delete adjacent_elements;adjacent_elements=0;delete faces;faces=0;
    delete node_on_boundary;node_on_boundary=0;delete boundary_mesh;boundary_mesh=0;delete face_hexahedrons;face_hexahedrons=0;
}
//#####################################################################
// Function Initialize_Incident_Elements
//#####################################################################
void HEXAHEDRON_MESH::
Initialize_Incident_Elements()
{
    delete incident_elements;incident_elements=new ARRAY<ARRAY<int> >(number_nodes);
    for(int h=1;h<=elements.m;h++){ // for each hexahedron, put it on each of its nodes lists of hexahedrons
        int p1,p2,p3,p4,p5,p6,p7,p8;elements(h).Get(p1,p2,p3,p4,p5,p6,p7,p8);
        (*incident_elements)(p1).Append(h);(*incident_elements)(p2).Append(h);
        (*incident_elements)(p3).Append(h);(*incident_elements)(p4).Append(h);
        (*incident_elements)(p5).Append(h);(*incident_elements)(p6).Append(h);
        (*incident_elements)(p7).Append(h);(*incident_elements)(p8).Append(h);}
    for(int i=1;i<=number_nodes;i++) (*incident_elements)(i).Compact(); // remove extra space
}
//#####################################################################
// Function Initialize_Adjacent_Elements
//#####################################################################
void HEXAHEDRON_MESH::
Initialize_Adjacent_Elements()
{
    delete adjacent_elements;adjacent_elements=new ARRAY<ARRAY<int> >(elements.m);
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int h=1;h<=elements.m;h++){ // for each hexahedron, make the list of adjacent hexahedrons
        int p1,p2,p3,p4,p5,p6,p7,p8;elements(h).Get(p1,p2,p3,p4,p5,p6,p7,p8);
        Find_And_Append_Adjacent_Elements(h,p1,p2,p3,p4);Find_And_Append_Adjacent_Elements(h,p5,p6,p7,p8);Find_And_Append_Adjacent_Elements(h,p3,p4,p7,p8);
        Find_And_Append_Adjacent_Elements(h,p1,p2,p5,p6);Find_And_Append_Adjacent_Elements(h,p2,p4,p6,p8);Find_And_Append_Adjacent_Elements(h,p1,p3,p5,p7);}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
    for(int h=1;h<=elements.m;h++) (*adjacent_elements)(h).Compact(); // remove extra space
}
//#####################################################################
// Function Find_And_Append_Adjacent_Elements
//#####################################################################
// check for an adjacent hexahedron that contains node1, node2, node3 and node4 - append to the adjaceny list
void HEXAHEDRON_MESH::
Find_And_Append_Adjacent_Elements(const int hexahedron,const int node1,const int node2,const int node3,const int node4)
{
    for(int h=1;h<=(*incident_elements)(node1).m;h++){
        int hexahedron2=(*incident_elements)(node1)(h);
        if(hexahedron2!=hexahedron && Node_In_Hexahedron(node2,hexahedron2) && Node_In_Hexahedron(node3,hexahedron2) && Node_In_Hexahedron(node4,hexahedron2))
            (*adjacent_elements)(hexahedron).Append_Unique(hexahedron2);}
}
//#####################################################################
// Function Initialize_Faces
//#####################################################################
void HEXAHEDRON_MESH::
Initialize_Faces()
{
    delete faces;faces=new ARRAY<VECTOR<int,4> >;
    HASHTABLE<VECTOR<int,4> > quad_list(2*elements.m); // list of faces currently found
    for(int h=1;h<=elements.m;h++){
        const VECTOR<int,8>& nodes=elements(h);
        for(int f=0;f<6;f++){
            VECTOR<int,4> face;for(int k=1;k<=4;k++) face[k]=nodes(face_indices[f][k-1]);
            if(quad_list.Set(face.Sorted())) faces->Append(VECTOR<int,4>(nodes(face[1]),nodes(face[2]),nodes(face[3]),nodes(face[4])));}}
}
//#####################################################################
// Function Initialize_Node_On_Boundary
//#####################################################################
void HEXAHEDRON_MESH::
Initialize_Node_On_Boundary()
{
    delete node_on_boundary;node_on_boundary=new ARRAY<bool>(number_nodes);ARRAY<int> p(8);
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int h=1;h<=elements.m;h++){
        elements(h).Get(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8));
        for(int f=0;f<6;f++)if(Number_Of_Hexahedrons_Across_Face(h,p(face_indices[f][0]),p(face_indices[f][1]),p(face_indices[f][2]),p(face_indices[f][3])) == 0){
            (*node_on_boundary)(p(face_indices[f][0]))=true;(*node_on_boundary)(p(face_indices[f][1]))=true;
            (*node_on_boundary)(p(face_indices[f][2]))=true;(*node_on_boundary)(p(face_indices[f][3]))=true;}}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Initialize_Boundary_Nodes
//#####################################################################
void HEXAHEDRON_MESH::
Initialize_Boundary_Nodes()
{
    delete boundary_nodes;boundary_nodes=new ARRAY<int>;ARRAY<int> p(8);
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int h=1;h<=elements.m;h++){
        elements(h).Get(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8));
        for(int f=0;f<6;f++)if(Number_Of_Hexahedrons_Across_Face(h,p(face_indices[f][0]),p(face_indices[f][1]),p(face_indices[f][2]),p(face_indices[f][3])) == 0){
            boundary_nodes->Append(p(face_indices[f][0]));boundary_nodes->Append(p(face_indices[f][1]));
            boundary_nodes->Append(p(face_indices[f][2]));boundary_nodes->Append(p(face_indices[f][3]));}}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}

    Sort(*boundary_nodes);
    boundary_nodes->Prune_Duplicates();
}
//#####################################################################
// Function Number_Of_Hexahedrons_Across_Face
//#####################################################################
int HEXAHEDRON_MESH::
Number_Of_Hexahedrons_Across_Face(const int hexahedron,const int node1,const int node2,const int node3,const int node4) const
{
    assert(incident_elements);int count=0;ARRAY<int> p(8);
    int n1=node1,n2=node2,n3=node3,n4=node4;exchange_sort(n1,n2,n3,n4);
    for(int h=1;h<=(*incident_elements)(n1).m;h++){
        int hexahedron2=(*incident_elements)(n1)(h);if(hexahedron==hexahedron2) continue; // hexahedron in question
        elements(hexahedron2).Get(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8));
        for(int f=0;f<6;f++){
            int i[4];for(int k=0;k<4;k++)i[k]=p(face_indices[f][k]);exchange_sort(i[0],i[1],i[2],i[3]);
            if(n1==i[0] && n2==i[1] && n3==i[2] && n4==i[3]){count++;break;}}}
    return count;
}
//#####################################################################
// Function Delete_Hexahedrons_With_Missing_Nodes
//#####################################################################
int HEXAHEDRON_MESH::
Delete_Hexahedrons_With_Missing_Nodes()
{
    ARRAY<int> deletion_list(elements.m);
    int number_deleted=0;
    for(int t=1;t<=elements.m;t++){
        int i[8];elements(t).Get(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]);
        if(!(i[0]&&i[1]&&i[2]&&i[3]&&i[4]&&i[5]&&i[6]&&i[7])) deletion_list(++number_deleted)=t;}
    deletion_list.Resize(number_deleted);
    Delete_Hexahedrons(deletion_list);
    return number_deleted;
}
//#####################################################################
// Function Delete_Hexahedrons
//#####################################################################
void HEXAHEDRON_MESH::
Delete_Hexahedrons(const ARRAY<int>& deletion_list)
{
    int last_index=elements.m;
    for(int k=1;k<=deletion_list.m;k++){
        int index=deletion_list(k);
        if(index > last_index) index=elements(index)(1); // hexahedron was moved, reset index to its new location
        for(int kk=1;kk<=8;kk++) elements(index)(kk)=elements(last_index)(kk);
        elements(last_index)(1)=index; // remembers where a hexahedron was moved to, in case it needs to be deleted later
        last_index--;}
    elements.Exact_Resize(last_index);
    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Initialize_Boudary_Mesh
//#####################################################################
void HEXAHEDRON_MESH::
Initialize_Boundary_Mesh()
{
    delete boundary_mesh;ARRAY<VECTOR<int,3> > triangle_list;
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    for(int h=1;h<=elements.m;h++)for(int f=0;f<6;f++){
        int p[4];for(int k=0;k<4;k++)p[k]=elements(h)(face_indices[f][k]);
        if(Number_Of_Hexahedrons_Across_Face(h,p[0],p[1],p[2],p[3]) == 0){
            triangle_list.Append(VECTOR<int,3>(p[0],p[1],p[2]));triangle_list.Append(VECTOR<int,3>(p[0],p[2],p[3]));}}
    boundary_mesh=new TRIANGLE_MESH(number_nodes,triangle_list);
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Initialize_Face_Hexahedrons
//#####################################################################
void HEXAHEDRON_MESH::
Initialize_Face_Hexahedrons()
{
    if(!faces) Initialize_Faces();
    bool incident_elements_defined=incident_elements!=0;if(!incident_elements) Initialize_Incident_Elements();
    delete face_hexahedrons;face_hexahedrons=new ARRAY<VECTOR<int,2> >(faces->m);
    for(int f=1;f<=faces->m;f++){
        int node1=(*faces)(f)(1),node2=(*faces)(f)(2),node3=(*faces)(f)(3),node4=(*faces)(f)(4),count=0;exchange_sort(node1,node2,node3,node4);
        for(int h=1;h<=(*incident_elements)(node1).m;h++){
            for(int hf=0;hf<6;hf++){
                ARRAY<int> p(4);for(int k=0;k<4;k++)p(k+1)=elements((*incident_elements)(node1)(h))(face_indices[hf][k]);exchange_sort(p(1),p(2),p(3),p(4));
                if(node1==p(1) && node2==p(2) && node3==p(3) && node4==p(4)){(*face_hexahedrons)(f)(++count)=(*incident_elements)(node1)(h);break;}}
            if(count>1) break;}}
    if(!incident_elements_defined){delete incident_elements;incident_elements=0;}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
void HEXAHEDRON_MESH::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int t=1;t<=elements.m;t++) for(int i=1;i<=7;i++) for(int j=i+1;j<=8;j++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(elements(t)[i],elements(t)[j]));
}
//#####################################################################
// Function Set_Number_Nodes
//#####################################################################
void HEXAHEDRON_MESH::
Set_Number_Nodes(const int number_nodes_input)
{
    if(number_nodes_input<number_nodes) PHYSBAM_FATAL_ERROR();
    number_nodes=number_nodes_input;
    if(incident_elements) incident_elements->Resize(number_nodes);
    if(node_on_boundary) node_on_boundary->Resize(number_nodes);
    if(boundary_mesh) boundary_mesh->Set_Number_Nodes(number_nodes);
}
//#####################################################################
// Function Mark_Nodes_Referenced
//#####################################################################
void HEXAHEDRON_MESH::
Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const
{
    for(int e=1;e<=elements.m;e++){INDIRECT_ARRAY<ARRAY<int>,VECTOR<int,8>&> subset=marks.Subset(elements(e));ARRAYS_COMPUTATIONS::Fill(subset,mark);}
}
//#####################################################################
