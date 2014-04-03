//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#ifndef __RIGID_GEOMETRY_RASTERIZATION_UNIFORM__
#define __RIGID_GEOMETRY_RASTERIZATION_UNIFORM__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>

namespace PhysBAM{
template<class TV> class COLLISION_GEOMETRY;
template<class T_GRID, class ID> class OBJECTS_IN_CELL;
template<class TV> class RANGE;

namespace RASTERIZATION{
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Object(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID& id);
template<class T,class TV,class T_GRID> void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,
    typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE& occupied,const bool with_body_motion,const T& extra_thickness,const T& body_thickness_factor);
template<class TV,class T_GRID> void Rasterize_Box_Onto_Blocks(const T_GRID& grid,typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE& occupied,const RANGE<TV>& box);
template<class TV,class T_GRID> void Rasterize_Box(const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const RANGE<TV>& box,const COLLISION_GEOMETRY_ID id);
//#####################################################################
};
};
#endif
