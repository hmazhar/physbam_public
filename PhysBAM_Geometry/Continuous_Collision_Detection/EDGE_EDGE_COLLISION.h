//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EDGE_EDGE_COLLISION__
#define __EDGE_EDGE_COLLISION__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T> class POINT_2D;
template<class T> class SEGMENT_3D;

namespace CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS
{
template<class T,class T_ARRAY> bool
Edge_Edge_Collision(const POINT_2D<T>& pt,const POINT_2D<T>& point,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V_edges,const T dt,const T collision_thickness,T& collision_time,
    VECTOR<T,2>& normal,VECTOR<T,2>& weights,T& relative_speed,bool allow_negative_weights,const T small_number,const bool exit_early);

template<class T> bool
Edge_Edge_Collision(const SEGMENT_3D<T>& seg_fault,const SEGMENT_3D<T>& segment,const VECTOR<T,3>& v1,const VECTOR<T,3>& v2,const VECTOR<T,3>& v3,const VECTOR<T,3>& v4,const T dt,
    const T collision_thickness,T& collision_time,VECTOR<T,3>& normal,VECTOR<T,2>& weights,T& relative_speed,const T small_number,const bool exit_early);
}
}
#endif
