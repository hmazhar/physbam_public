//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TETRAHEDRALIZED_VOLUME_CENTROID__
#define __TETRAHEDRALIZED_VOLUME_CENTROID__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
namespace PhysBAM{
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Centroid_Of_Neighbors
//#####################################################################
template<class T>
VECTOR<T,3> Centroid_Of_Neighbors(const TETRAHEDRALIZED_VOLUME<T>& tv,const int node)
{
    typedef VECTOR<T,3> TV;
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool neighbor_nodes_defined=tv.mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined) tv.mesh.Initialize_Neighbor_Nodes();
    TV target;
    int number_of_neighbors=(*tv.mesh.neighbor_nodes)(node).m;
    for(int j=1;j<=number_of_neighbors;j++) target+=tv.particles.X((*tv.mesh.neighbor_nodes)(node)(j));
    if(number_of_neighbors != 0) target/=(T)number_of_neighbors;
    else target=tv.particles.X(node); // if no neighbors, return the current node location
    if(!neighbor_nodes_defined){delete tv.mesh.neighbor_nodes;tv.mesh.neighbor_nodes=0;}
    return target;
}
}
}
#endif
