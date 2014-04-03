//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TETRAHEDRALIZED_VOLUME_SPLIT__
#define __TETRAHEDRALIZED_VOLUME_SPLIT__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
namespace PhysBAM{
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Split_Along_Fracture_Plane
//#####################################################################
template<class T>
void Split_Along_Fracture_Plane(TETRAHEDRALIZED_VOLUME<T>& tv,const PLANE<T>& plane,ARRAY<int>& particle_replicated)
{
    typedef VECTOR<T,3> TV;
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    if(particle_replicated.m != tv.particles.array_collection->Size()) particle_replicated.Resize(tv.particles.array_collection->Size());
    bool incident_elements_defined=tv.mesh.incident_elements!=0;if(!incident_elements_defined) tv.mesh.Initialize_Incident_Elements();
    ARRAY<bool> positive_side(tv.mesh.elements.m);int number_on_positive_side=0;
    int t;for(t=1;t<=tv.mesh.elements.m;t++){int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
        TV x1=tv.particles.X(i),x2=tv.particles.X(j),x3=tv.particles.X(k),x4=tv.particles.X(l),centroid=(T).25*(x1+x2+x3+x4);
        if(plane.Signed_Distance(centroid) >= 0){positive_side(t)=true;number_on_positive_side++;}}
    int p;for(p=1;p<=tv.particles.array_collection->Size();p++){
        bool seen_positive_side=false,seen_negative_side=false;
        for(t=1;t<=(*tv.mesh.incident_elements)(p).m;t++){
            if(positive_side((*tv.mesh.incident_elements)(p)(t))) seen_positive_side=true;else seen_negative_side=true;}
        if(seen_positive_side && seen_negative_side) particle_replicated(p)=1;}
    int number_of_new_particles=0;
    for(p=1;p<=tv.particles.array_collection->Size();p++) if(particle_replicated(p)){ // assumes we're storing mass (this is not set here!), position, & velocity
        int new_index=tv.particles.array_collection->Add_Element();particle_replicated(p)=new_index;number_of_new_particles++;
        tv.particles.X(new_index)=tv.particles.X(p);tv.particles.V(new_index)=tv.particles.V(p);}
    // loop through tets and change indices in negative_side to new_indices (from particle_replicated)
    for(t=1;t<=tv.mesh.elements.m;t++) if(!positive_side(t)){ int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l); // replace indices with replicated_indices
        if(particle_replicated(i)) i=particle_replicated(i);if(particle_replicated(j)) j=particle_replicated(j);
        if(particle_replicated(k)) k=particle_replicated(k);if(particle_replicated(l)) l=particle_replicated(l);
        tv.mesh.elements(t).Set(i,j,k,l);}
    tv.mesh.number_nodes=tv.particles.array_collection->Size();
    if(incident_elements_defined){delete tv.mesh.incident_elements;tv.mesh.incident_elements=0;}
}
//#####################################################################
// Function Split_Node
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements
template<class T>
int Split_Node(TETRAHEDRALIZED_VOLUME<T>& tv,const int particle_index,const VECTOR<T,3>& normal)
{
    typedef VECTOR<T,3> TV;
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool incident_elements_defined=tv.mesh.incident_elements!=0;if(!incident_elements_defined) tv.mesh.Initialize_Incident_Elements();
    PLANE<T> plane(normal,tv.particles.X(particle_index));ARRAY<int> tets_incident_on_old_particle,tets_incident_on_new_particle;
    int t;for(t=1;t<=(*tv.mesh.incident_elements)(particle_index).m;t++){
        int this_incident_tet=(*tv.mesh.incident_elements)(particle_index)(t);int i,j,k,l;tv.mesh.elements(this_incident_tet).Get(i,j,k,l);
        TV x1=tv.particles.X(i),x2=tv.particles.X(j),x3=tv.particles.X(k),x4=tv.particles.X(l),centroid=(T).25*(x1+x2+x3+x4);
        if(plane.Signed_Distance(centroid) < 0) tets_incident_on_new_particle.Append(this_incident_tet);
        else tets_incident_on_old_particle.Append(this_incident_tet);}
    int new_particle=0;
    if(tets_incident_on_old_particle.m != 0 && tets_incident_on_new_particle.m != 0){
        // new particle - assumes we're storing position, and velocity - user must fix mass outside this function call
        new_particle=tv.particles.array_collection->Add_Element();tv.mesh.number_nodes=tv.particles.array_collection->Size();
        tv.particles.X(new_particle)=tv.particles.X(particle_index);tv.particles.V(new_particle)=tv.particles.V(particle_index);
        for(t=1;t<=(*tv.mesh.incident_elements)(particle_index).m;t++){
            int this_incident_tet=(*tv.mesh.incident_elements)(particle_index)(t);int i,j,k,l;tv.mesh.elements(this_incident_tet).Get(i,j,k,l);
            TV x1=tv.particles.X(i),x2=tv.particles.X(j),x3=tv.particles.X(k),x4=tv.particles.X(l),centroid=(T).25*(x1+x2+x3+x4);
            if(plane.Signed_Distance(centroid) < 0){ // relabel with duplicate node
                if(i == particle_index) i=new_particle;if(j == particle_index) j=new_particle;if(k == particle_index) k=new_particle;if(l == particle_index) l=new_particle;
                tv.mesh.elements(this_incident_tet).Set(i,j,k,l);}}
        if(incident_elements_defined){ //repair incident tetrahedrons if necessary
            (*tv.mesh.incident_elements)(particle_index).Clean_Memory();
            (*tv.mesh.incident_elements)(particle_index).Append_Elements(tets_incident_on_old_particle);
            (*tv.mesh.incident_elements).Append(tets_incident_on_new_particle);}}
    if(!incident_elements_defined){delete tv.mesh.incident_elements;tv.mesh.incident_elements=0;}
    return new_particle;
}
//#####################################################################
// Function Split_Connected_Component
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements & adjacent_elements
template<class T>
int Split_Connected_Component(TETRAHEDRALIZED_VOLUME<T>& tv,const int node)
{
    typedef VECTOR<T,3> TV;
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool incident_elements_defined=tv.mesh.incident_elements!=0;if(!incident_elements_defined) tv.mesh.Initialize_Incident_Elements();
    bool adjacent_elements_defined=tv.mesh.adjacent_elements!=0;if(!adjacent_elements_defined) tv.mesh.Initialize_Adjacent_Elements();
    ARRAY<bool> marked((*tv.mesh.incident_elements)(node).m);
    tv.mesh.Mark_Face_Connected_Component_Incident_On_A_Node(node,marked);
    int number_marked=0;int t;for(t=1;t<=marked.m;t++) if(marked(t)) number_marked++;
    int new_particle=0;
    if(number_marked != marked.m){
        // new particle -- assumes we're storing position, and velocity - user must fix mass outside this function call
        new_particle=tv.particles.array_collection->Add_Element();tv.mesh.number_nodes=tv.particles.array_collection->Size();
        tv.particles.X(new_particle)=tv.particles.X(node);tv.particles.V(new_particle)=tv.particles.V(node);
        ARRAY<int> empty;tv.mesh.incident_elements->Append(empty);
        ARRAY<int> indices_to_remove(number_marked);int counter=0;
        for(t=1;t<=(*tv.mesh.incident_elements)(node).m;t++) if(marked(t)){
            indices_to_remove(++counter)=t;
            for(int i=1;i<=4;i++) if(tv.mesh.elements((*tv.mesh.incident_elements)(node)(t))(i) == node){
                tv.mesh.elements((*tv.mesh.incident_elements)(node)(t))(i)=new_particle;break;}
            (*tv.mesh.incident_elements)(new_particle).Append((*tv.mesh.incident_elements)(node)(t));}
        (*tv.mesh.incident_elements)(node).Remove_Sorted_Indices(indices_to_remove);}
    if(!incident_elements_defined){delete tv.mesh.incident_elements;tv.mesh.incident_elements=0;}
    if(!adjacent_elements_defined){delete tv.mesh.adjacent_elements;tv.mesh.adjacent_elements=0;}
    return new_particle;
}
}
}
#endif
