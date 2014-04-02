//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGLE_SUBDIVISION_LINEAR__
#define __TRIANGLE_SUBDIVISION_LINEAR__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Apply_Linear_Subdivision
//#####################################################################
template<class TV>
void Apply_Linear_Subdivision(TRIANGLE_SUBDIVISION& ts,ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==ts.start_index_for_new_nodes-1+ts.triangle_mesh.segment_mesh->elements.m);
    // neighbor_nodes are used to identify which nodes are in the base mesh - construct here and delete when this class is deleted
    if(!ts.triangle_mesh.neighbor_nodes){ts.delete_neighbor_nodes=true;ts.triangle_mesh.Initialize_Neighbor_Nodes();}
    for(int i=1;i<=ts.triangle_mesh.number_nodes;i++)if((*ts.triangle_mesh.neighbor_nodes)(i).m) subdivided_values(i)=base_values(i);
    // interpolate values on edges
    for(int k=1;k<=ts.triangle_mesh.segment_mesh->elements.m;k++)
        subdivided_values(ts.start_index_for_new_nodes-1+k)=(T).5*(base_values(ts.triangle_mesh.segment_mesh->elements(k)(1))+base_values(ts.triangle_mesh.segment_mesh->elements(k)(2)));
}
}
}
#endif
