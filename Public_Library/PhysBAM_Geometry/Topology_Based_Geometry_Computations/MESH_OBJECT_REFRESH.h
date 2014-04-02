//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MESH_OBJECT_REFRESH__
#define __MESH_OBJECT_REFRESH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>

namespace PhysBAM{
template<class TV,class T_MESH> class MESH_OBJECT;
template<class TV> class BOX;
template<class TV> class RANGE;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
template<class TV,class T_MESH> void
Update_Bounding_Box(MESH_OBJECT<TV,T_MESH>& mo)
{
    if(!mo.bounding_box) mo.bounding_box=new BOX<TV>();
    const BOX_HIERARCHY<TV>* hierarchy=&*debug_cast<const typename MESH_TO_OBJECT<TV,T_MESH>::TYPE&>(mo).hierarchy;
    if(!mo.mesh.elements.m) *mo.bounding_box=RANGE<TV>::Bounding_Box(mo.particles.X);
    else if(hierarchy) *mo.bounding_box=hierarchy->box_hierarchy(hierarchy->root);
    else *mo.bounding_box=RANGE<TV>::Bounding_Box(mo.particles.X.Subset(mo.mesh.elements.Flattened()));
}
//#####################################################################
// Function Initialize_Particle_Partition
//#####################################################################
template<class TV,class T_MESH> void
Initialize_Particle_Partition(MESH_OBJECT<TV,T_MESH>& mo,const VECTOR<int,TV::m>& counts)
{
    PHYSBAM_ASSERT(mo.bounding_box);VECTOR<int,TV::m> counts_new;
    for(int i=1;i<=counts.m;i++) counts_new[i]=mo.desired_particle_partition_counts[i]?mo.desired_particle_partition_counts[i]:counts[i];
    PHYSBAM_ASSERT(counts_new.All_Greater(VECTOR<int,TV::m>()));
    delete mo.particle_partition;mo.particle_partition=new PARTICLE_PARTITION<TV>(*mo.bounding_box,counts_new,mo.particles);
}
}
}
#endif
