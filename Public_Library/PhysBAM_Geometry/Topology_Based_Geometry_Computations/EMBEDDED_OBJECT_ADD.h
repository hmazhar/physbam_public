//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EMBEDDED_OBJECT_ADD__
#define __EMBEDDED_OBJECT_ADD__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_OBJECT.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Add_Embedded_Particle_To_Embedded_Children
//#####################################################################
template<class TV,int d> void
Add_Embedded_Particle_To_Embedded_Children(EMBEDDED_OBJECT<TV,d>& eo,const int embedded_particle)
{
    int parent1,parent2;eo.parent_particles(embedded_particle).Get(parent1,parent2);
    int index1=(*eo.embedded_children_index)(parent1),index2=(*eo.embedded_children_index)(parent2);
    if(index1) (*eo.embedded_children)(index1).Append(embedded_particle);
    else{
        eo.embedded_children->Append(ARRAY<int>());(*eo.embedded_children_index)(parent1)=eo.embedded_children->m;
        (*eo.embedded_children)(eo.embedded_children->m).Append(embedded_particle);}
    if(index2) (*eo.embedded_children)(index2).Append(embedded_particle);
    else{
        eo.embedded_children->Append(ARRAY<int>());(*eo.embedded_children_index)(parent2)=eo.embedded_children->m;
        (*eo.embedded_children)(eo.embedded_children->m).Append(embedded_particle);}
}
//#####################################################################
// Function Add_Emdedded_Particle
//#####################################################################
template<class TV,class T,int d> int
Add_Embedded_Particle(EMBEDDED_OBJECT<TV,d>& eo,const VECTOR<int,2>& nodes,const T interpolation_fraction_input,const bool reinitialize_hash_table)
{
    assert(eo.simplicial_object.mesh.number_nodes==eo.particles.array_collection->Size() && eo.embedded_mesh.number_nodes==eo.particles.array_collection->Size() && !eo.Embedded_Particle_On_Segment(nodes));
    int new_embedded_particle=eo.embedded_particles.Add_Element();
    eo.simplicial_object.mesh.Add_Nodes(1);eo.embedded_mesh.Add_Nodes(1); // Need to update meshes and acceleration structures
    eo.interpolation_fraction.Append(eo.Clamp_Interpolation_Fraction(interpolation_fraction_input));
    eo.parent_particles.Append(nodes);
    eo.parents_to_embedded_particles_hash_table->Insert(nodes.Sorted(),new_embedded_particle);
    if(reinitialize_hash_table) eo.Initialize_Parents_To_Embedded_Particles_Hash_Table();
    if(eo.embedded_children_index) eo.Add_Embedded_Particle_To_Embedded_Children(new_embedded_particle);
    T lambda=eo.interpolation_fraction(new_embedded_particle);
    eo.embedded_particles.X(new_embedded_particle)=LINEAR_INTERPOLATION<T,TV>::Linear(eo.particles.X(nodes[1]),eo.particles.X(nodes[2]),lambda);
    return new_embedded_particle;
}
//#####################################################################
// Function Add_Emdedded_Particle_If_Not_Already_There
//#####################################################################
template<class TV,class T,int d> int
Add_Embedded_Particle_If_Not_Already_There(EMBEDDED_OBJECT<TV,d>& eo,const VECTOR<int,2>& nodes,const T interpolation_fraction_input)
{
    if(int current_embedded_particle=eo.Embedded_Particle_On_Segment(nodes)) return current_embedded_particle;
    else return eo.Add_Embedded_Particle(nodes,interpolation_fraction_input);
}
//#####################################################################
// Function Add_Embedded_Subelement_To_Embedded_Subelements_In_Element
//#####################################################################
template<class TV,int d> void
Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(EMBEDDED_OBJECT<TV,d>& eo,const int subelement)
{
    int element=eo.Element_Containing_Subelement(subelement);
    int index=(*eo.embedded_subelements_in_parent_element_index)(element);
    if(index){
        int embedded_subelement_number=++(*eo.number_of_embedded_subelements_in_parent_element)(index);
        (*eo.embedded_subelements_in_parent_element)(index)(embedded_subelement_number)=subelement;}
    else{
        VECTOR<int,EMBEDDED_OBJECT<TV,d>::max_subelements_per_element> subelements;subelements[1]=subelement;
        eo.embedded_subelements_in_parent_element->Append(subelements);
        eo.number_of_embedded_subelements_in_parent_element->Append(1);
        (*eo.embedded_subelements_in_parent_element_index)(element)=eo.embedded_subelements_in_parent_element->m;}
}
//#####################################################################
// Function Add_Embedded_Subelement_If_Not_Already_There
//#####################################################################
template<class TV,int d> int
Add_Embedded_Subelement_If_Not_Already_There(EMBEDDED_OBJECT<TV,d>& eo,const VECTOR<int,d>& embedded_nodes)
{
    VECTOR<int,d> global_particles(eo.embedded_particles.active_indices.Subset(embedded_nodes));
    if(int t=eo.embedded_mesh.Simplex(global_particles)) return t;
    else return eo.Add_Embedded_Subelement(embedded_nodes);
}
//#####################################################################
// Function Add_Embedded_Subelement
//#####################################################################
// only updates embedded_mesh.incident_elements and embedded_subelements_in_parent_element
template<class TV,int d> int
Add_Embedded_Subelement(EMBEDDED_OBJECT<TV,d>& eo,const VECTOR<int,d>& embedded_nodes)
{
    VECTOR<int,d> global_particles(eo.embedded_particles.active_indices.Subset(embedded_nodes));
    assert(!eo.embedded_mesh.Simplex(global_particles));
    int new_subelement=eo.embedded_mesh.elements.Append(global_particles);
    if(eo.embedded_mesh.incident_elements) // needs to be updated
        for(int i=1;i<=global_particles.m;i++) (*eo.embedded_mesh.incident_elements)(global_particles[i]).Append(new_subelement);
    if(eo.embedded_subelements_in_parent_element) eo.Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(new_subelement);
    return new_subelement;
}
}
}
#endif
