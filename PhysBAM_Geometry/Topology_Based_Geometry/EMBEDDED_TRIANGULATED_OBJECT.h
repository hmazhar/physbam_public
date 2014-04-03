//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TRIANGULATED_OBJECT
//#####################################################################
#ifndef __EMBEDDED_TRIANGULATED_OBJECT__
#define __EMBEDDED_TRIANGULATED_OBJECT__

#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV>
class EMBEDDED_TRIANGULATED_OBJECT:public EMBEDDED_OBJECT<TV,2>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE T_SEGMENTED_CURVE;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
public:
    typedef EMBEDDED_OBJECT<TV,2> BASE;
    using BASE::embedded_particles;using BASE::parent_particles;using BASE::Is_Parent;using BASE::Embedded_Particle_On_Segment;using BASE::embedded_mesh;using BASE::simplicial_object;
    using BASE::Embedded_Subelements_In_Element;using BASE::Number_Of_Embedded_Subelements_In_Element;using BASE::Node_In_Simplex_Is_Material;using BASE::Add_Embedded_Subelement;

    EMBEDDED_TRIANGULATED_OBJECT(T_TRIANGULATED_OBJECT& simplicial_object_input);

    int Node_Separated_By_Embedded_Subelement(const int embedded_segment) const
    {int global_a,global_b;embedded_mesh.elements(embedded_segment).Get(global_a,global_b);
    int a=embedded_particles.subset_index_from_point_cloud_index(global_a),b=embedded_particles.subset_index_from_point_cloud_index(global_b);
    int ppa1,ppa2,ppb1,ppb2;parent_particles(a).Get(ppa1,ppa2);parent_particles(b).Get(ppb1,ppb2);
    if(Is_Parent(ppa1,b))return ppa1;else return ppa2;}

    bool Nodes_Are_Separated_In_Simplex(const int node1,const int node2,const int triangle) const
    {int embedded_node=Embedded_Particle_On_Segment(node1,node2);if(!embedded_node) return false;
    VECTOR<int,2> segments=Embedded_Subelements_In_Element(triangle);
    int global_embedded_node=embedded_particles.active_indices(embedded_node);
    if(!segments[1]) return false;else if(embedded_mesh.Node_In_Segment(global_embedded_node,segments[1])) return true;
    if(!segments[2]) return false;else if(embedded_mesh.Node_In_Segment(global_embedded_node,segments[2])) return true;
    return false;}

    bool Nodes_Are_Materially_Connected_In_Simplex(const int node1,const int node2,const int simplex) const
    {return Node_In_Simplex_Is_Material(node1,simplex) && Node_In_Simplex_Is_Material(node2,simplex) && !Nodes_Are_Separated_In_Simplex(node1,node2,simplex);}

    int Number_Of_Embedded_Cuts(const int triangle) // a quad only counts as one cut
    {return Number_Of_Embedded_Subelements_In_Element(triangle);}
   
    bool Segment_Is_Broken(const int node1,const int node2) const
    {if(!Embedded_Particle_On_Segment(node1,node2)) return false;
    ARRAY<int> triangles_on_edge;simplicial_object.mesh.Triangles_On_Edge(node1,node2,&triangles_on_edge);
    for(int t=1;t<=triangles_on_edge.m;t++) if(!Nodes_Are_Separated_In_Simplex(node1,node2,triangles_on_edge(t)))return false;
    return true;}

    int Number_Of_Edges_With_Embedded_Particles(const int triangle)
    {int i,j,k;simplicial_object.mesh.elements(triangle).Get(i,j,k);
    int ij=Embedded_Particle_On_Segment(i,j),ik=Embedded_Particle_On_Segment(i,k),jk=Embedded_Particle_On_Segment(j,k);
    return (ij>0)+(ik>0)+(jk>0);}

    int Embedded_Node_Common_To_Both_Segments_In_Triangle(const int triangle)
    {VECTOR<int,2> emb_segments=Embedded_Subelements_In_Element(triangle);assert(emb_segments[1] && emb_segments[2]);
    int global_a,global_b,global_c,global_d;embedded_mesh.elements(emb_segments[1]).Get(global_a,global_b);embedded_mesh.elements(emb_segments[2]).Get(global_c,global_d);
    return embedded_particles.subset_index_from_point_cloud_index((global_a==global_c || global_a==global_d)?global_a:global_b);}

    int Isolated_Node(const int triangle)
    {VECTOR<int,2> emb_segments=Embedded_Subelements_In_Element(triangle);assert(emb_segments[1] && !emb_segments[2]);
    return Node_Separated_By_Embedded_Subelement(emb_segments[1]);}
    
    int Diamond_Node(const int triangle)
    {if(Number_Of_Embedded_Subelements_In_Element(triangle)<2) return 0;
    int common_emb_node=Embedded_Node_Common_To_Both_Segments_In_Triangle(triangle);
    int i,j,k;simplicial_object.mesh.elements(triangle).Get(i,j,k);
    if(!Is_Parent(i,common_emb_node))return i;
    if(!Is_Parent(j,common_emb_node))return j;
    return k;}

    T Fraction_Of_Triangles_With_N_Cuts(const int n)
    {int count=0;for(int t=1;t<=simplicial_object.mesh.elements.m;t++) if(Number_Of_Embedded_Cuts(t) == n) count++;
    return (T)count/(T)simplicial_object.mesh.elements.m;}

    int Add_Embedded_Segment(const int embedded_particle1,const int embedded_particle2)
    {return Add_Embedded_Subelement(VECTOR<int,2>(embedded_particle1,embedded_particle2));}

    int Add_Embedded_Segment(const int embedded_particle1_parent1,const int embedded_particle1_parent2,const int embedded_particle2_parent1,
        const int embedded_particle2_parent2)
    {return Add_Embedded_Segment(Embedded_Particle_On_Segment(embedded_particle1_parent1,embedded_particle1_parent2),
        Embedded_Particle_On_Segment(embedded_particle2_parent1,embedded_particle2_parent2));}

//#####################################################################
};
}
#endif
