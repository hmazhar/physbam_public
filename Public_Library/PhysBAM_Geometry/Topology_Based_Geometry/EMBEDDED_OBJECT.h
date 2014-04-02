//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_OBJECT
//#####################################################################
#ifndef __EMBEDDED_OBJECT__
#define __EMBEDDED_OBJECT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_SUBSET.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

static VECTOR<int,2> endpoints_temp;

template<class TV,int d>
class EMBEDDED_OBJECT:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d-1>::OBJECT T_BOUNDARY_OBJECT;
    typedef typename MESH_POLICY<d-1>::MESH T_BOUNDARY_MESH;
    typedef typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
public:
    enum WORKAROUND1 {topological_dimension=d};
    enum WORKAROUND2 {max_subelements_per_element=2*d-2};

    GEOMETRY_PARTICLES<TV>& particles; // reference to the particles array containing all particles (including embedded)
    POINT_CLOUD_SUBSET<TV,GEOMETRY_PARTICLES<TV> > embedded_particles;
    ARRAY<VECTOR<int,2> > parent_particles; // has length embedded_particles.active_indices.m, indexes into particles
    ARRAY<T> interpolation_fraction; // has length embedded_particles.active_indices.m
    T interpolation_fraction_threshold;
    ARRAY<int>* embedded_children_index; // has length particles.array_collection->Size(), indexes into embedded children
    ARRAY<ARRAY<int> >* embedded_children; // indexes into embedded_particles
    HASHTABLE<VECTOR<int,2>,int>* parents_to_embedded_particles_hash_table; // indexes into embedded_particles
    int hashtable_multiplier;
    ARRAY<int>* embedded_subelements_in_parent_element_index; // length number of parent elements
    ARRAY<int>* number_of_embedded_subelements_in_parent_element; // lower dimensional 
    ARRAY<VECTOR<int,max_subelements_per_element> >* embedded_subelements_in_parent_element;
    bool average_interpolation_fractions;

    T_SIMPLICIAL_OBJECT& simplicial_object;
    T_BOUNDARY_MESH embedded_mesh; // indexes into particles
    T_BOUNDARY_OBJECT embedded_object;
    ARRAY<VECTOR<bool,d+1> > node_in_simplex_is_material;
    ARRAY<int> orientation_index; // currently only used sometimes in 3d
    BOX<TV>* bounding_box;
protected:
    bool need_destroy_simplicial_object;
public:

protected:
    EMBEDDED_OBJECT(T_SIMPLICIAL_OBJECT& simplicial_object_input);
public:
    virtual ~EMBEDDED_OBJECT();

    void Set_Interpolation_Fraction_Threshold(const T interpolation_fraction_threshold_input=.1)
    {interpolation_fraction_threshold=interpolation_fraction_threshold_input;}

    T Clamp_Interpolation_Fraction(const T interpolation_fraction_input) const
    {return clamp(interpolation_fraction_input,interpolation_fraction_threshold,(T)1-interpolation_fraction_threshold);} 
    
    int Embedded_Particle_On_Segment(const int endpoint1,const int endpoint2) const
    {return Embedded_Particle_On_Segment(VECTOR<int,2>(endpoint1,endpoint2));}

    int Embedded_Particle_On_Segment(const VECTOR<int,2>& endpoints) const
    {endpoints_temp=endpoints; // compiler bug, gcc 4.0.1
    int embedded_particle=0;parents_to_embedded_particles_hash_table->Get(endpoints.Sorted(),embedded_particle);
    return embedded_particle;}

    int Particle_Embedded_On_Segment(const int endpoint1,const int endpoint2) const
    {int embedded_particle=Embedded_Particle_On_Segment(endpoint1,endpoint2);
    return embedded_particle?embedded_particles.active_indices(embedded_particle):0;}

    bool Is_Parent(const int parent_node,const int embedded_node) const
    {return parent_particles(embedded_node).Contains(parent_node);}

    int Which_Parent(const int parent_node,const int embedded_node) const
    {return parent_particles(embedded_node).Find(parent_node);}

    bool Are_Parents(const VECTOR<int,2>& parents,const int embedded_node) const
    {return Is_Parent(parents[1],embedded_node) && Is_Parent(parents[2],embedded_node);}

    int Other_Parent(const int parent,const int embedded_node) const
    {int index=Which_Parent(parent,embedded_node);if(index == 1) return parent_particles(embedded_node)(2);else if(index == 2) return parent_particles(embedded_node)(1);else return 0;}

    void Replace_Parent_Particle(const int embedded_particle,const int old_parent_particle,const int new_parent_particle)
    {if(parent_particles(embedded_particle)(1) == old_parent_particle) parent_particles(embedded_particle)(1)=new_parent_particle;
    else{assert(parent_particles(embedded_particle)(2) == old_parent_particle);parent_particles(embedded_particle)(2)=new_parent_particle;}}

    int Number_Of_Children(const int parent_node) const
    {if(!(*embedded_children_index)(parent_node)) return 0;else return (*embedded_children)((*embedded_children_index)(parent_node)).m;}

    int Child(const int parent_node,const int child_index) const
    {return (*embedded_children)((*embedded_children_index)(parent_node))(child_index);}

    T Fraction_Of_Elements_With_Embedded_Subelements()
    {if(!embedded_subelements_in_parent_element_index) return 0;
    int count=0;for(int t=1;t<=embedded_subelements_in_parent_element_index->m;t++) if((*embedded_subelements_in_parent_element_index)(t)) count++;
    return (T)count/(T)embedded_subelements_in_parent_element_index->m;}

    void Update_Bounding_Box()
    {simplicial_object.Update_Bounding_Box();*bounding_box=*simplicial_object.bounding_box;}

    void Initialize_Node_In_Simplex_Is_Material()
    {if(node_in_simplex_is_material.m==simplicial_object.mesh.elements.m) return; // don't overwrite existing values if they exist
    node_in_simplex_is_material.Resize(simplicial_object.mesh.elements.m,false,false);
    ARRAY_VIEW<bool> view=node_in_simplex_is_material.Flattened();ARRAYS_COMPUTATIONS::Fill(view,true);}

    bool Node_In_Simplex_Is_Material(const int node,const int element) const
    {int index=simplicial_object.mesh.elements(element).Find(node);assert(index);
    return node_in_simplex_is_material(element)(index);}

    bool Node_Near_Material(const int node) const
    {const ARRAY<int>& incident=(*simplicial_object.mesh.incident_elements)(node);
    for(int i=1;i<=incident.m;i++) if(Node_In_Simplex_Is_Material(node,incident(i))) return true;
    return false;}

    void Initialize_Orientation_Index_If_Necessary()
    {if(!orientation_index.m) orientation_index.Resize(simplicial_object.mesh.elements.m);
    else if(orientation_index.m!=simplicial_object.mesh.elements.m) PHYSBAM_FATAL_ERROR();} // already initialized to the wrong size

    VECTOR<int,max_subelements_per_element> Embedded_Subelements_In_Element(const int element) const
    {assert(embedded_subelements_in_parent_element);
    int index=(*embedded_subelements_in_parent_element_index)(element);
    return index?(*embedded_subelements_in_parent_element)(index):VECTOR<int,max_subelements_per_element>();}

    int Number_Of_Embedded_Subelements_In_Element(const int element) const
    {assert(embedded_subelements_in_parent_element);
    int index=(*embedded_subelements_in_parent_element_index)(element);
    return index?(*number_of_embedded_subelements_in_parent_element)(index):0;}

    TV Position_Of_Embedded_Particle(const int embedded_particle) const
    {int pp1,pp2;parent_particles(embedded_particle).Get(pp1,pp2);
    return LINEAR_INTERPOLATION<T,TV>::Linear(particles.X(pp1),particles.X(pp2),interpolation_fraction(embedded_particle));}

    TV Position_Of_Embedded_Particle(const int i,const int j) const
    {return Position_Of_Embedded_Particle(Embedded_Particle_On_Segment(i,j));}

    int Add_Embedded_Particle_If_Not_Already_There(const int node1,const int node2,const T interpolation_fraction_input)
    {return Add_Embedded_Particle_If_Not_Already_There(VECTOR<int,2>(node1,node2),interpolation_fraction_input);}

    void Update_Number_Nodes() PHYSBAM_OVERRIDE
    {simplicial_object.Update_Number_Nodes();embedded_mesh.Set_Number_Nodes(particles.array_collection->Size());embedded_particles.Update_Number_Nodes();}

    void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const PHYSBAM_OVERRIDE
    {simplicial_object.Mark_Nodes_Referenced(marks,mark);embedded_object.Mark_Nodes_Referenced(marks,mark);}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("EMBEDDED_OBJECT<T,VECTOR<T,%d>,%d>",TV::dimension,d);}

//#####################################################################
    virtual void Clean_Memory();
    void Delete_Auxiliary_Structures();
    static T_EMBEDDED_OBJECT* Create();
    static T_EMBEDDED_OBJECT* Create(GEOMETRY_PARTICLES<TV>& particles);
    void Reset_Parents_To_Embedded_Particles_Hash_Table();
    void Copy_Then_Reset_Parents_To_Embedded_Particles_Hash_Table(HASHTABLE<VECTOR<int,2>,int>*& hash_table_copy);
    void Initialize_Parents_To_Embedded_Particles_Hash_Table(const int new_hashtable_multiplier=0);
    void Update_Embedded_Particle_Positions();
    void Initialize_Embedded_Children();
    void Add_Embedded_Particle_To_Embedded_Children(const int embedded_particle);
    int Add_Embedded_Particle_If_Not_Already_There(const VECTOR<int,2>& nodes,const T interpolation_fraction_input);
    int Add_Embedded_Particle(const VECTOR<int,2>& nodes,const T interpolation_fraction_input,const bool reinitialize_hash_table=true);
    int Element_Containing_Subelement(const int embedded_subelement) const;
    void Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(const int subelement); // incrementally update embedded_subelements_in_parent_element
    int Add_Embedded_Subelement_If_Not_Already_There(const VECTOR<int,d>& embedded_particles);
    int Add_Embedded_Subelement(const VECTOR<int,d>& embedded_particle);
    void Initialize_Embedded_Subelements_In_Parent_Element();
    void Calculate_Boundary_From_Levelset_On_Nodes(ARRAY<T>& phi,const bool discard_elements_outside_levelset=true,const bool verbose=true);
//#####################################################################
};
}
#endif
