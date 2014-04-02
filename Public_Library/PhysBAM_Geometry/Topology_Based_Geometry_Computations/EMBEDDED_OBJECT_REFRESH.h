//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EMBEDDED_OBJECT_REFRESH__
#define __EMBEDDED_OBJECT_REFRESH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_OBJECT.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Initialize_Parents_To_Embedded_Particles_Hash_Table
//#####################################################################
template<class TV,int d> void
Initialize_Parents_To_Embedded_Particles_Hash_Table(EMBEDDED_OBJECT<TV,d>& eo,const int new_hashtable_multiplier)
{
    if(new_hashtable_multiplier) eo.hashtable_multiplier=new_hashtable_multiplier;
    if(eo.parents_to_embedded_particles_hash_table && eo.parents_to_embedded_particles_hash_table->Max_Size() >= eo.hashtable_multiplier*eo.particles.array_collection->Size()) return;
    delete eo.parents_to_embedded_particles_hash_table;eo.parents_to_embedded_particles_hash_table=new HASHTABLE<VECTOR<int,2>,int>(eo.hashtable_multiplier*eo.particles.array_collection->Size());
    for(int p=1;p<=eo.embedded_particles.active_indices.m;p++)
        eo.parents_to_embedded_particles_hash_table->Insert(eo.parent_particles(p).Sorted(),p);
}
//#####################################################################
// Function Update_Embedded_Particle_Positions
//#####################################################################
template<class TV,int d> void
Update_Embedded_Particle_Positions(EMBEDDED_OBJECT<TV,d>& eo)
{
    for(int p=1;p<=eo.embedded_particles.active_indices.m;p++){
        int i,j;eo.parent_particles(p).Get(i,j);
        eo.embedded_particles.X(p)=LINEAR_INTERPOLATION<typename TV::SCALAR,TV>::Linear(eo.particles.X(i),eo.particles.X(j),eo.interpolation_fraction(p));}
}
//#####################################################################
// Function Initialize_Embedded_Children
//#####################################################################
template<class TV,int d> void
Initialize_Embedded_Children(EMBEDDED_OBJECT<TV,d>& eo)
{
    delete eo.embedded_children_index;delete eo.embedded_children;
    eo.embedded_children_index=new ARRAY<int>(eo.particles.array_collection->Size());eo.embedded_children=new ARRAY<ARRAY<int> >();
    for(int p=1;p<=eo.embedded_particles.active_indices.m;p++) eo.Add_Embedded_Particle_To_Embedded_Children(p);
}
//#####################################################################
// Function Initialize_Embedded_Subelements_In_Parent_Element
//#####################################################################
template<class TV,int d> void
Initialize_Embedded_Subelements_In_Parent_Element(EMBEDDED_OBJECT<TV,d>& eo)
{
    delete eo.embedded_subelements_in_parent_element_index;delete eo.number_of_embedded_subelements_in_parent_element;delete eo.embedded_subelements_in_parent_element;
    eo.embedded_subelements_in_parent_element_index=new ARRAY<int>(eo.simplicial_object.mesh.elements.m);
    eo.number_of_embedded_subelements_in_parent_element=new ARRAY<int>(0);
    eo.embedded_subelements_in_parent_element=new ARRAY<VECTOR<int,EMBEDDED_OBJECT<TV,d>::max_subelements_per_element> >();
    for(int t=1;t<=eo.embedded_mesh.elements.m;t++) eo.Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(t);
}
template<class TV> void Add_Levelset_Cuts(EMBEDDED_TRIANGULATED_OBJECT<TV>& embedded_object,const int positive_count,const VECTOR<int,3>& nodes)
{
    int i,j,k;nodes.Get(i,j,k);
    if(positive_count==1) embedded_object.Add_Embedded_Segment(i,k,j,k); // one segment
    else if(positive_count==2) embedded_object.Add_Embedded_Segment(i,j,i,k);
}
template<class T> void Add_Levelset_Cuts(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object,const int positive_count,const VECTOR<int,4>& nodes)
{
    int i,j,k,l;nodes.Get(i,j,k,l);
    if(positive_count==1) embedded_object.Add_Embedded_Triangle(i,l,j,l,k,l); // one triangle
    else if(positive_count==2){ // two triangles
        int ik=embedded_object.Embedded_Particle_On_Segment(i,k),il=embedded_object.Embedded_Particle_On_Segment(i,l),
            jk=embedded_object.Embedded_Particle_On_Segment(j,k),jl=embedded_object.Embedded_Particle_On_Segment(j,l);
        embedded_object.Add_Embedded_Triangle(il,jl,jk);embedded_object.Add_Embedded_Triangle(jk,ik,il);}
    else if(positive_count==3) embedded_object.Add_Embedded_Triangle(i,j,i,k,i,l); // one triangle
}
//#####################################################################
// Function Calculate_Boundary_From_Levelset_On_Nodes
//#####################################################################
template<class TV,class T,int d> void
Calculate_Boundary_From_Levelset_On_Nodes(EMBEDDED_OBJECT<TV,d>& eo,ARRAY<T>& phi,const bool discard_elements_outside_levelset,const bool verbose)
{
    if(discard_elements_outside_levelset){ // TODO: Consider doing this before initializing the embedded object
        assert(eo.embedded_particles.active_indices.m==0); // Discarding should be done before adding any embedded particles
        for(int t=eo.simplicial_object.mesh.elements.m;t>0;t--){
            VECTOR<T,d+1> phi_t(phi.Subset(eo.simplicial_object.mesh.elements(t)));
            if(phi_t.Min()>0) eo.simplicial_object.mesh.elements.Remove_Index_Lazy(t);}
        ARRAY<int> condensation_mapping;eo.simplicial_object.Discard_Valence_Zero_Particles_And_Renumber(condensation_mapping);
        ARRAY<T> phi_new(eo.particles.array_collection->Size());for(int k=1;k<=phi.m;k++) if(condensation_mapping(k)) phi_new(condensation_mapping(k))=phi(k);
        phi.Exchange(phi_new);}

    if(eo.embedded_particles.active_indices.m){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Calculate_Boundary_From_Levelset_On_Nodes cannot be called with existing embedded particles"<<std::endl;
#endif
        PHYSBAM_FATAL_ERROR();}

    eo.Initialize_Parents_To_Embedded_Particles_Hash_Table();
    eo.embedded_particles.Update_Subset_Index_From_Element_Index(); // TODO: move this somewhere else
    eo.embedded_mesh.number_nodes=eo.particles.array_collection->Size();

    bool segment_mesh_defined=eo.simplicial_object.mesh.segment_mesh!=0;if(!segment_mesh_defined) eo.simplicial_object.mesh.Initialize_Segment_Mesh();
    bool embedded_incident_elements_defined=eo.embedded_mesh.incident_elements!=0;if(!embedded_incident_elements_defined) eo.embedded_mesh.Initialize_Incident_Elements();

    // calculate embedded particles
    for(int s=1;s<=eo.simplicial_object.mesh.segment_mesh->elements.m;s++){
        int n1,n2;eo.simplicial_object.mesh.segment_mesh->elements(s).Get(n1,n2);T phi1=phi(n1),phi2=phi(n2);
        if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
            int inside,outside;if(phi1 <= 0){inside=n1;outside=n2;}else{inside=n2;outside=n1;exchange(phi1,phi2);}
            T theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
            eo.Add_Embedded_Particle(VECTOR<int,2>(inside,outside),theta);}}

    // calculate embedded simplices
    for(int t=1;t<=eo.simplicial_object.mesh.elements.m;t++){
        VECTOR<int,d+1> nodes=eo.simplicial_object.mesh.elements(t);
        {int i=1,j=nodes.m;while(i<j){if(phi(nodes[i])>0) exchange(nodes[i],nodes[j--]);else i++;} // move inside nodes before outside nodes
        if((nodes.m-j)&1){if(phi(nodes[2])<=0) exchange(nodes[1],nodes[2]);else exchange(nodes[d],nodes[d+1]);}} // one final swap to ensure an even permutation
        for(int i=1;i<nodes.m;i++) assert((phi(nodes[i])>0) <= (phi(nodes[i+1])>0));
        int positive_count=0;for(int i=1;i<=nodes.m;i++) if(phi(nodes[i])>0) positive_count++;
        Add_Levelset_Cuts(dynamic_cast<typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT&>(eo),positive_count,nodes);}

    eo.node_in_simplex_is_material.Resize(eo.simplicial_object.mesh.elements.m);
    for(int t=1;t<=eo.simplicial_object.mesh.elements.m;t++){VECTOR<int,d+1>& element=eo.simplicial_object.mesh.elements(t);
        for(int i=1;i<=element.m;i++) eo.node_in_simplex_is_material(t)(i)=phi(element[i])<=0;}

    if(!segment_mesh_defined){delete eo.simplicial_object.mesh.segment_mesh;eo.simplicial_object.mesh.segment_mesh=0;}
    if(!embedded_incident_elements_defined){delete eo.embedded_mesh.incident_elements;eo.embedded_mesh.incident_elements=0;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(verbose){
        LOG::cout << "total particles = " << eo.particles.array_collection->Size() << std::endl;
        LOG::cout << "total elements = " << eo.simplicial_object.mesh.elements.m << std::endl;
        LOG::cout << "total embedded subelements = " << eo.embedded_mesh.elements.m << std::endl;}
#endif
}
}
}
#endif
