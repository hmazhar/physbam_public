//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TETRAHEDRALIZED_VOLUME_PRUNE__
#define __TETRAHEDRALIZED_VOLUME_PRUNE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
namespace PhysBAM{
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Discard_Spikes_From_Adjacent_Elements
//#####################################################################
// throws out all tetrahedrons with only one neighbor (i.e. a spike on the boundary)
// returns index of first discarded
template<class T>
void Discard_Spikes_From_Adjacent_Elements(TETRAHEDRALIZED_VOLUME<T>& tv,ARRAY<int>* deletion_list)
{
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=tv.mesh.adjacent_elements!=0;if(!adjacent_elements_defined) tv.mesh.Initialize_Adjacent_Elements();
    if(deletion_list){
        deletion_list->Resize(0);
        for(int t=1;t<=tv.mesh.elements.m;t++) if((*tv.mesh.adjacent_elements)(t).m == 1) deletion_list->Append(t);
        tv.mesh.Delete_Sorted_Elements(*deletion_list);}
    else{
        ARRAY<int> list;
        for(int t=1;t<=tv.mesh.elements.m;t++) if((*tv.mesh.adjacent_elements)(t).m == 1) list.Append(t);
        tv.mesh.Delete_Sorted_Elements(list);}
    if(!adjacent_elements_defined){delete tv.mesh.adjacent_elements;tv.mesh.adjacent_elements=0;}
}
//#####################################################################
// Function Interior_Edges_With_Boundary_Nodes
//#####################################################################
// throws out all tetrahedrons with all four nodes on boundary
template<class T>
void Interior_Edges_With_Boundary_Nodes(TETRAHEDRALIZED_VOLUME<T>& tv,ARRAY<VECTOR<int,2> >* deletion_list)
{
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool neighbor_nodes_defined=tv.mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined) tv.mesh.Initialize_Neighbor_Nodes();
    bool node_on_boundary_defined=tv.mesh.node_on_boundary!=0;if(!node_on_boundary_defined) tv.mesh.Initialize_Node_On_Boundary();
    bool boundary_nodes_defined=tv.mesh.boundary_nodes!=0;if(!boundary_nodes_defined) tv.mesh.Initialize_Boundary_Nodes();
    bool boundary_mesh_defined=tv.mesh.boundary_mesh!=0;if(!boundary_mesh_defined) tv.mesh.Initialize_Boundary_Mesh();
    bool boundary_mesh_segment_mesh_defined=tv.mesh.boundary_mesh->segment_mesh!=0;
    if(!boundary_mesh_segment_mesh_defined) tv.mesh.boundary_mesh->Initialize_Segment_Mesh();
    bool segment_mesh_defined=tv.mesh.segment_mesh!=0;if(!segment_mesh_defined) tv.mesh.Initialize_Segment_Mesh();
    bool boundary_mesh_segment_mesh_incident_elements_defined=tv.mesh.boundary_mesh->segment_mesh->incident_elements!=0;
    if(!boundary_mesh_segment_mesh_incident_elements_defined) tv.mesh.boundary_mesh->segment_mesh->Initialize_Incident_Elements();

    deletion_list->Resize(0);
    for(int t=1;t<=tv.mesh.boundary_nodes->m;t++){
        int node1=(*tv.mesh.boundary_nodes)(t);
        for(int i=1;i<=(*tv.mesh.neighbor_nodes)(node1).m;i++){
            int node2=(*tv.mesh.neighbor_nodes)(node1)(i);
            if(node1 <node2 && (*tv.mesh.node_on_boundary)(node2) && !tv.mesh.boundary_mesh->segment_mesh->Segment(node1,node2))
                deletion_list->Append(VECTOR<int,2>(node1,node2));}}

    // delete node_on_boundary if defined in this function
    if(!boundary_mesh_segment_mesh_incident_elements_defined){
        delete tv.mesh.boundary_mesh->segment_mesh->incident_elements;tv.mesh.boundary_mesh->segment_mesh->incident_elements=0;}
    if(!segment_mesh_defined){delete tv.mesh.segment_mesh;tv.mesh.segment_mesh=0;}
    if(!boundary_mesh_segment_mesh_defined){delete tv.mesh.boundary_mesh->segment_mesh;tv.mesh.boundary_mesh->segment_mesh=0;}
    if(!boundary_mesh_defined){delete tv.mesh.boundary_mesh;tv.mesh.boundary_mesh=0;}
    if(!boundary_nodes_defined){delete tv.mesh.boundary_nodes;tv.mesh.boundary_nodes=0;}
    if(!node_on_boundary_defined){delete tv.mesh.node_on_boundary;tv.mesh.node_on_boundary=0;}
    if(!neighbor_nodes_defined){delete tv.mesh.neighbor_nodes;tv.mesh.neighbor_nodes=0;}
}
//#####################################################################
// Function Discard_Spikes
//#####################################################################
// throws out all tetrahedrons with all four nodes on boundary
template<class T>
void Discard_Spikes(TETRAHEDRALIZED_VOLUME<T>& tv,ARRAY<int>* deletion_list)
{
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool node_on_boundary_defined=tv.mesh.node_on_boundary!=0;if(!node_on_boundary_defined) tv.mesh.Initialize_Node_On_Boundary();
    if(deletion_list){
        deletion_list->Resize(0);
        for(int t=1;t<=tv.mesh.elements.m;t++){int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
            if((*tv.mesh.node_on_boundary)(i) && (*tv.mesh.node_on_boundary)(j) && (*tv.mesh.node_on_boundary)(k) && (*tv.mesh.node_on_boundary)(l))
                deletion_list->Append(t);}
        tv.mesh.Delete_Sorted_Elements(*deletion_list);}
    else{
        ARRAY<int> list;
        for(int t=1;t<=tv.mesh.elements.m;t++){int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
            if((*tv.mesh.node_on_boundary)(i) && (*tv.mesh.node_on_boundary)(j) && (*tv.mesh.node_on_boundary)(k) && (*tv.mesh.node_on_boundary)(l))
                list.Append(t);}
        tv.mesh.Delete_Sorted_Elements(list);}
    if(!node_on_boundary_defined){delete tv.mesh.node_on_boundary;tv.mesh.node_on_boundary=0;}
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface
//#####################################################################
// uses Whitney-like criterion to discard only those tets that are for sure outside levelset (assuming accurate signed distance)
template<class T>
void Discard_Tetrahedrons_Outside_Implicit_Surface(TETRAHEDRALIZED_VOLUME<T>& tv,IMPLICIT_OBJECT<VECTOR<T,3> >& implicit_surface)
{
    typedef VECTOR<T,3> TV;
    for(int t=tv.mesh.elements.m;t>=1;t--){
        int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
        TV xi=tv.particles.X(i),xj=tv.particles.X(j),xk=tv.particles.X(k),xl=tv.particles.X(l);
        T max_length=TETRAHEDRON<T>::Maximum_Edge_Length(xi,xj,xk,xl);
        T min_phi=min(implicit_surface(xi),implicit_surface(xj),implicit_surface(xk),implicit_surface(xl));
        if(min_phi>max_length) tv.mesh.elements.Remove_Index_Lazy(t);}
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive
//#####################################################################
// uses Whitney-like criterion to discard only those tets that are not for sure inside levelset (assuming accurate signed distance)
template<class T>
void Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(TETRAHEDRALIZED_VOLUME<T>& tv,IMPLICIT_OBJECT<VECTOR<T,3> >& implicit_surface)
{
    typedef VECTOR<T,3> TV;
    for(int t=tv.mesh.elements.m;t>=1;t--){
        int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
        TV xi=tv.particles.X(i),xj=tv.particles.X(j),xk=tv.particles.X(k),xl=tv.particles.X(l);
        T max_length=TETRAHEDRON<T>::Maximum_Edge_Length(xi,xj,xk,xl);
        T max_phi=max(implicit_surface(xi),implicit_surface(xj),implicit_surface(xk),implicit_surface(xl));
        if(max_phi+max_length>0) tv.mesh.elements.Remove_Index_Lazy(t);}
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive
//#####################################################################
// same as above Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive but restricted to tets fully contained inside one of the supplied boxes
// TODO: Need to decide what constitutes good criteria for tet to be contained in a box: for now demanding full containment, i.e. all four particles
template<class T>
void Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(TETRAHEDRALIZED_VOLUME<T>& tv,IMPLICIT_OBJECT<VECTOR<T,3> >& implicit_surface,const ARRAY<RANGE<VECTOR<T,3> > >& bounding_boxes)
{
    typedef VECTOR<T,3> TV;
    for(int b=1;b<=bounding_boxes.Size();b++){const RANGE<TV>& box=bounding_boxes(b);
        for(int t=tv.mesh.elements.m;t>=1;t--){
            int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
            TV xi=tv.particles.X(i),xj=tv.particles.X(j),xk=tv.particles.X(k),xl=tv.particles.X(l);
            if(box.Lazy_Outside(xi)||box.Lazy_Outside(xj)||box.Lazy_Outside(xk)||box.Lazy_Outside(xl)) continue;
            T max_length=TETRAHEDRON<T>::Maximum_Edge_Length(xi,xj,xk,xl);
            T max_phi=max(implicit_surface(xi),implicit_surface(xj),implicit_surface(xk),implicit_surface(xl));
            if(max_phi+max_length>0) tv.mesh.elements.Remove_Index_Lazy(t);}}
}
}
}
#endif
