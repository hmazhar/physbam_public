//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __POINT_FACE_COLLISION__
#define __POINT_FACE_COLLISION__

#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T> class SEGMENT_2D;
template<class T> class TRIANGLE_3D;

namespace CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS
{
template<class T> POINT_SIMPLEX_COLLISION_TYPE
Robust_Point_Triangle_Collision(const TRIANGLE_3D<T>& initial_triangle,const TRIANGLE_3D<T>& final_triangle,const VECTOR<T,3>& x,const VECTOR<T,3>& final_x,const T dt,
    const T collision_thickness,T& collision_time,VECTOR<T,3>& normal,VECTOR<T,3>& weights,T& relative_speed);
template<class T> bool
Point_Face_Collision(const TRIANGLE_3D<T>& tri,const VECTOR<T,3>& x,const VECTOR<T,3>& v,const VECTOR<T,3>& v1,const VECTOR<T,3>& v2,const VECTOR<T,3>& v3,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,3>& normal,VECTOR<T,3>& weights,T& relative_speed,const bool exit_early);

template<class T> POINT_SIMPLEX_COLLISION_TYPE
Robust_Point_Segment_Collision(const SEGMENT_2D<T>& initial_segment,const SEGMENT_2D<T>& final_segment,const VECTOR<T,2> &x,const VECTOR<T,2> &final_x,const T dt, const T collision_thickness,
    T& collision_time,VECTOR<T,2>& normal,T& collision_alpha,T& relative_speed);
template<class T> bool
Point_Face_Collision(const SEGMENT_2D<T>& seg_fault,const VECTOR<T,2>& x,const VECTOR<T,2>& v,const VECTOR<T,2>& v1,const VECTOR<T,2>& v2,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,2>& normal,VECTOR<T,2>& weights,T& relative_speed,const bool exit_early);

template<class T> POINT_SIMPLEX_COLLISION_TYPE
Robust_Point_Point_Collision(const POINT_SIMPLEX_1D<T>& initial_simplex,const POINT_SIMPLEX_1D<T>& final_simplex,const VECTOR<T,1>& x,const VECTOR<T,1>& final_x,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,1>& normal,ONE& weights,T& relative_speed);
template<class T> bool
Point_Point_Collision(const POINT_SIMPLEX_1D<T>& initial_simplex,const VECTOR<T,1>& x,const VECTOR<T,1>& v,const VECTOR<T,1>& v1,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,1>& normal,ONE& weights,T& relative_speed,const bool exit_early);
}
}
#endif
