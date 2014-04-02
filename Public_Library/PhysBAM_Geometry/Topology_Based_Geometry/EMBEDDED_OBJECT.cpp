//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_OBJECT
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/permutation.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/EMBEDDED_OBJECT_ADD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/EMBEDDED_OBJECT_REFRESH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> EMBEDDED_OBJECT<TV,d>::
EMBEDDED_OBJECT(T_SIMPLICIAL_OBJECT& simplicial_object_input)
    :particles(simplicial_object_input.particles),embedded_particles(particles),interpolation_fraction_threshold((T).1),
    embedded_children_index(0),embedded_children(0),parents_to_embedded_particles_hash_table(0),hashtable_multiplier(12),
    embedded_subelements_in_parent_element_index(0),number_of_embedded_subelements_in_parent_element(0),
    embedded_subelements_in_parent_element(0),average_interpolation_fractions(false),
    simplicial_object(simplicial_object_input),embedded_object(embedded_mesh,particles),bounding_box(0),need_destroy_simplicial_object(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> EMBEDDED_OBJECT<TV,d>::
~EMBEDDED_OBJECT()
{
    if(need_destroy_simplicial_object) delete &simplicial_object;
    delete embedded_children_index;delete embedded_children;delete parents_to_embedded_particles_hash_table;delete embedded_subelements_in_parent_element_index;
    delete number_of_embedded_subelements_in_parent_element;delete embedded_subelements_in_parent_element;delete bounding_box;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Clean_Memory()
{
    embedded_particles.active_indices.Clean_Memory();parent_particles.Clean_Memory();interpolation_fraction.Clean_Memory();
    embedded_mesh.Clean_Memory();node_in_simplex_is_material.Clean_Memory();Delete_Auxiliary_Structures();
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Delete_Auxiliary_Structures()
{
    delete embedded_children_index;embedded_children_index=0;delete embedded_children;embedded_children=0;
    delete parents_to_embedded_particles_hash_table;parents_to_embedded_particles_hash_table=0;
    delete embedded_subelements_in_parent_element_index;embedded_subelements_in_parent_element_index=0;
    delete number_of_embedded_subelements_in_parent_element;number_of_embedded_subelements_in_parent_element=0;
    delete embedded_subelements_in_parent_element;embedded_subelements_in_parent_element=0;
    embedded_object.Clean_Memory();delete bounding_box;bounding_box=0;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT* EMBEDDED_OBJECT<TV,d>::
Create()
{
    T_EMBEDDED_OBJECT* object=new T_EMBEDDED_OBJECT(*T_SIMPLICIAL_OBJECT::Create());
    object->need_destroy_simplicial_object=true;return object;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT* EMBEDDED_OBJECT<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    T_EMBEDDED_OBJECT* object=new T_EMBEDDED_OBJECT(*T_SIMPLICIAL_OBJECT::Create(particles));
    object->Initialize_Parents_To_Embedded_Particles_Hash_Table(15); // TODO: reconsider
    object->Set_Interpolation_Fraction_Threshold((T)1e-4); // TODO: reconsider
    object->need_destroy_simplicial_object=true;return object;
}
//#####################################################################
// Function Reset_Parents_To_Embedded_Particles_Hash_Table
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Reset_Parents_To_Embedded_Particles_Hash_Table()
{
    if(parents_to_embedded_particles_hash_table) parents_to_embedded_particles_hash_table->Remove_All();
}
//#####################################################################
// Function Copy_Then_Reset_Parents_To_Embedded_Particles_Hash_Table
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Copy_Then_Reset_Parents_To_Embedded_Particles_Hash_Table(HASHTABLE<VECTOR<int,2>,int>*& hash_table_copy)
{
    hash_table_copy=parents_to_embedded_particles_hash_table;parents_to_embedded_particles_hash_table=new HASHTABLE<VECTOR<int,2>,int>(hashtable_multiplier*particles.array_collection->Size()); // TODO: Rethink the size
}
//#####################################################################
// Function Initialize_Parents_To_Embedded_Particles_Hash_Table
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Initialize_Parents_To_Embedded_Particles_Hash_Table(const int new_hashtable_multiplier)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Parents_To_Embedded_Particles_Hash_Table(*this,new_hashtable_multiplier);
}
//#####################################################################
// Function Update_Embedded_Particle_Positions
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Update_Embedded_Particle_Positions()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Embedded_Particle_Positions(*this);
}
//#####################################################################
// Function Initialize_Embedded_Children
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Initialize_Embedded_Children()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Embedded_Children(*this);
}
//#####################################################################
// Function Add_Embedded_Particle_To_Embedded_Children
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Particle_To_Embedded_Children(const int embedded_particle)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Add_Embedded_Particle_To_Embedded_Children(*this,embedded_particle);
}
//#####################################################################
// Function Add_Emdedded_Particle_If_Not_Already_There
//#####################################################################
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Particle_If_Not_Already_There(const VECTOR<int,2>& nodes,const T interpolation_fraction_input)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Add_Embedded_Particle_If_Not_Already_There(*this,nodes,interpolation_fraction_input);
}
//#####################################################################
// Function Add_Emdedded_Particle
//#####################################################################
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Particle(const VECTOR<int,2>& nodes,const T interpolation_fraction_input,const bool reinitialize_hash_table)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Add_Embedded_Particle(*this,nodes,interpolation_fraction_input,reinitialize_hash_table);
}
//#####################################################################
// Function Element_Containing_Subelement
//#####################################################################
template<int d> static inline VECTOR<int,d+1> Compact_Parents(const VECTOR<VECTOR<int,2>,d>& parents)
{
    VECTOR<int,d+1> all_parents(parents[1]);int m=2;
    for(int i=2;i<=d;i++)for(int j=1;j<=2;j++)
        if(!all_parents.Contains(parents[i][j])) all_parents[++m]=parents[i][j];
    assert(m==d+1);
    return all_parents;
}
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Element_Containing_Subelement(const int embedded_subelement) const
{
    VECTOR<int,d> embedded_nodes(embedded_particles.subset_index_from_point_cloud_index.Subset(embedded_mesh.elements(embedded_subelement)));
    VECTOR<VECTOR<int,2>,d> parents(parent_particles.Subset(embedded_nodes));
    return simplicial_object.mesh.Simplex(Compact_Parents(parents));
}
//#####################################################################
// Function Add_Embedded_Subelement_To_Embedded_Subelements_In_Element
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(const int subelement)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(*this,subelement);
}
//#####################################################################
// Function Add_Embedded_Subelement_If_Not_Already_There
//#####################################################################
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Subelement_If_Not_Already_There(const VECTOR<int,d>& embedded_nodes)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Add_Embedded_Subelement_If_Not_Already_There(*this,embedded_nodes);
}
//#####################################################################
// Function Add_Embedded_Subelement
//#####################################################################
// only updates embedded_mesh.incident_elements and embedded_subelements_in_parent_element
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Subelement(const VECTOR<int,d>& embedded_nodes)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Add_Embedded_Subelement(*this,embedded_nodes);
}
//#####################################################################
// Function Initialize_Embedded_Subelements_In_Parent_Element
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Initialize_Embedded_Subelements_In_Parent_Element()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Embedded_Subelements_In_Parent_Element(*this);
}
//#####################################################################
// Function Calculate_Boundary_From_Levelset_On_Nodes
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Calculate_Boundary_From_Levelset_On_Nodes(ARRAY<T>& phi,const bool discard_elements_outside_levelset,const bool verbose)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Calculate_Boundary_From_Levelset_On_Nodes(*this,phi,discard_elements_outside_levelset,verbose);
}
//#####################################################################
template class EMBEDDED_OBJECT<VECTOR<float,2>,2>;
template class EMBEDDED_OBJECT<VECTOR<float,3>,2>;
template class EMBEDDED_OBJECT<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EMBEDDED_OBJECT<VECTOR<double,2>,2>;
template class EMBEDDED_OBJECT<VECTOR<double,3>,2>;
template class EMBEDDED_OBJECT<VECTOR<double,3>,3>;
#endif
