//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header COLLISIONS_GEOMETRY_FORWARD
//#####################################################################
#ifndef __COLLISIONS_GEOMETRY_FORWARD__
#define __COLLISIONS_GEOMETRY_FORWARD__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class COLLISION_GEOMETRY;
template<class TV> class COLLISION_GEOMETRY_COLLECTION;
template<class TV> class COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR;

template<class T_COLLISION_GEOMETRY,class T_ARRAY,class ID=int> class COLLISION_GEOMETRY_SPATIAL_PARTITION;
enum SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC {SPATIAL_PARTITION_SCENE_SIZE,SPATIAL_PARTITION_MAX_BOX_SIZE,SPATIAL_PARTITION_AVERAGE_BOX_SIZE};
}
#endif
