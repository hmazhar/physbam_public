//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const SEGMENT_3D<T>& segment,const TRIANGLE_3D<T>& triangle,const T thickness_over_two)
{
    RAY<VECTOR<T,3> > ray(segment.x1,segment.x2-segment.x1);ray.semi_infinite=false;ray.t_max=(segment.x2-segment.x1).Magnitude();
    return INTERSECTION::Intersects(ray,triangle,thickness_over_two);
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const SEGMENT_3D<T>& segment,const TRIANGLE_3D<T>& triangle,T& a,VECTOR<T,3>& weights,const T thickness_over_two)
{
    RAY<VECTOR<T,3> > ray(segment.x1,segment.x2-segment.x1);ray.semi_infinite=false;
    T magnitude=(segment.x2-segment.x1).Magnitude();ray.t_max=magnitude;
    if(!INTERSECTION::Intersects(ray,triangle,thickness_over_two)){a=ray.t_max/magnitude;weights=triangle.Barycentric_Coordinates(ray.Point(ray.t_max));return true;}
    return false;
}
//#####################################################################
template bool Intersects(const SEGMENT_3D<float>&,const TRIANGLE_3D<float>&,const float);
template bool Intersects(const SEGMENT_3D<float>&,const TRIANGLE_3D<float>&,float&,VECTOR<float,3>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const SEGMENT_3D<double>&,const TRIANGLE_3D<double>&,const double);
template bool Intersects(const SEGMENT_3D<double>&,const TRIANGLE_3D<double>&,double&,VECTOR<double,3>&,const double);
#endif
};
};
