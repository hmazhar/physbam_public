//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_AREA_SPLIT__
#define __TRIANGULATED_AREA_SPLIT__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Split_Node
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements
template<class T>
int Split_Node(TRIANGULATED_AREA<T>& ta,const int particle_index,const VECTOR<T,2>& normal)
{
    int incident_elements_defined=ta.mesh.incident_elements!=0;if(!incident_elements_defined) ta.mesh.Initialize_Incident_Elements();
    VECTOR<T,2> x0=ta.particles.X(particle_index);ARRAY<int> tris_incident_on_old_particle,tris_incident_on_new_particle;
    int t;for(t=1;t<=(*ta.mesh.incident_elements)(particle_index).m;t++){
        int this_incident_tri=(*ta.mesh.incident_elements)(particle_index)(t);int i,j,k;ta.mesh.elements(this_incident_tri).Get(i,j,k);
        VECTOR<T,2> x1=ta.particles.X(i),x2=ta.particles.X(j),x3=ta.particles.X(k),centroid=(T)one_third*(x1+x2+x3);
        if((centroid.x-x0.x)*normal.x+(centroid.y-x0.y)*normal.y < 0) tris_incident_on_new_particle.Append(this_incident_tri); 
        else tris_incident_on_old_particle.Append(this_incident_tri);}
    int new_particle=0;
    if(tris_incident_on_old_particle.m != 0 && tris_incident_on_new_particle.m != 0){ 
        // new particle - assumes we're storing position, and velocity - user must fix mass outside this function call
        new_particle=ta.particles.array_collection->Add_Element();ta.mesh.number_nodes=ta.particles.array_collection->Size();
        ta.particles.X(new_particle)=ta.particles.X(particle_index);ta.particles.V(new_particle)=ta.particles.V(particle_index);
        for(t=1;t<=(*ta.mesh.incident_elements)(particle_index).m;t++){
            int this_incident_tri=(*ta.mesh.incident_elements)(particle_index)(t);int i,j,k;ta.mesh.elements(this_incident_tri).Get(i,j,k);
            VECTOR<T,2> x1=ta.particles.X(i),x2=ta.particles.X(j),x3=ta.particles.X(k),centroid=(T)one_third*(x1+x2+x3);
            if((centroid.x-x0.x)*normal.x+(centroid.y-x0.y)*normal.y < 0){ // relabel with duplicate node
                if(i == particle_index) i=new_particle;if(j == particle_index) j=new_particle;if(k == particle_index) k=new_particle;
                ta.mesh.elements(this_incident_tri).Set(i,j,k);}}        
        if(incident_elements_defined){ //repair incident triangles if necessary
            (*ta.mesh.incident_elements)(particle_index).Clean_Memory();
            (*ta.mesh.incident_elements)(particle_index).Append_Elements(tris_incident_on_old_particle);
            (*ta.mesh.incident_elements).Append(tris_incident_on_new_particle);}}
    if(!incident_elements_defined){delete ta.mesh.incident_elements;ta.mesh.incident_elements=0;}
    return new_particle;
}
//#####################################################################
// Function Split_Connected_Component
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements & adjacent_elements
template<class T>
int Split_Connected_Component(TRIANGULATED_AREA<T>& ta,const int node)
{
    int incident_elements_defined=ta.mesh.incident_elements!=0;if(!incident_elements_defined) ta.mesh.Initialize_Incident_Elements();
    int adjacent_elements_defined=ta.mesh.adjacent_elements!=0;if(!adjacent_elements_defined) ta.mesh.Initialize_Adjacent_Elements();
    ARRAY<bool> marked((*ta.mesh.incident_elements)(node).m);
    ta.mesh.Mark_Edge_Connected_Component_Incident_On_A_Node(node,marked);
    int number_marked=0;int t;for(t=1;t<=marked.m;t++) if(marked(t)) number_marked++;
    int new_particle=0;
    if(number_marked != marked.m){
        // new particle -- assumes we're storing position, and velocity - user must fix mass outside this function call
        new_particle=ta.particles.array_collection->Add_Element();ta.mesh.number_nodes=ta.particles.array_collection->Size();
        ta.particles.X(new_particle)=ta.particles.X(node);ta.particles.V(new_particle)=ta.particles.V(node);
        ARRAY<int> empty;ta.mesh.incident_elements->Append(empty);
        ARRAY<int> indices_to_remove(number_marked);int counter=0;
        for(t=1;t<=(*ta.mesh.incident_elements)(node).m;t++) if(marked(t)){
            indices_to_remove(++counter)=t;
            for(int i=1;i<=3;i++) if(ta.mesh.elements((*ta.mesh.incident_elements)(node)(t))(i) == node){ 
                ta.mesh.elements((*ta.mesh.incident_elements)(node)(t))(i)=new_particle;break;}
            (*ta.mesh.incident_elements)(new_particle).Append((*ta.mesh.incident_elements)(node)(t));}
        (*ta.mesh.incident_elements)(node).Remove_Sorted_Indices(indices_to_remove);}
    if(!incident_elements_defined){delete ta.mesh.incident_elements;ta.mesh.incident_elements=0;} 
    if(!adjacent_elements_defined){delete ta.mesh.adjacent_elements;ta.mesh.adjacent_elements=0;}
    return new_particle;
}
}
}
#endif
