//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGLE_SUBDIVISION_FRACTAL__
#define __TRIANGLE_SUBDIVISION_FRACTAL__

#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Apply_Fractal_Subdivision
//#####################################################################
template<class TV>
void Apply_Fractal_Subdivision(TRIANGLE_SUBDIVISION& ts,ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values,const float power)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==ts.start_index_for_new_nodes-1+ts.triangle_mesh.segment_mesh->elements.m);
    // neighbor_nodes are used to identify which nodes are in the base mesh - construct here and delete when this class is deleted
    if(!ts.triangle_mesh.neighbor_nodes){ts.delete_neighbor_nodes=true;ts.triangle_mesh.Initialize_Neighbor_Nodes();}
    if(!ts.triangle_mesh.incident_elements){ts.delete_incident_elements=true;ts.triangle_mesh.Initialize_Incident_Elements();}
    for(int i=1;i<=ts.triangle_mesh.number_nodes;i++)if((*ts.triangle_mesh.neighbor_nodes)(i).m) subdivided_values(i)=base_values(i);
    // interpolate values on edges
    RANDOM_NUMBERS<T> random;
    for(int k=1;k<=ts.triangle_mesh.segment_mesh->elements.m;k++){
        int node1=ts.triangle_mesh.segment_mesh->elements(k)(1);
        int node2=ts.triangle_mesh.segment_mesh->elements(k)(2);
        TV midpoint=(T).5*(base_values(node1)+base_values(node2));
        TV offset=base_values(node1)-base_values(node2);
        T modulus=(T)power*random.Get_Gaussian();
        ARRAY<int> adjacent_triangles;
        ts.triangle_mesh.Triangles_On_Edge(node1,node2,&adjacent_triangles);
        TV normal;
        for(int i=1;i<=adjacent_triangles.m;i++){
            int n1,n2,n3;ts.triangle_mesh.elements(adjacent_triangles(i)).Get(n1,n2,n3);
            normal+=TRIANGLE_3D<T>::Normal(base_values(n1),base_values(n2),base_values(n3));}
        normal.Normalize();
        //subdivided_values(start_index_for_new_nodes-1+k)=midpoint+normal*edge_magnitude*modulus;
        subdivided_values(ts.start_index_for_new_nodes-1+k)=midpoint+normal*modulus;
    }
}
}
}
#endif
