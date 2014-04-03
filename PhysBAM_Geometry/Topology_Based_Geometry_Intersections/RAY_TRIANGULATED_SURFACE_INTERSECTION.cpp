//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const TRIANGULATED_SURFACE<T>& surface, const T thickness_over_two)
{
    PHYSBAM_ASSERT(surface.hierarchy);
    return surface.hierarchy->Intersection(ray,thickness_over_two);
}
//#####################################################################
// Function Closest_Non_Intersecting_Point
//#####################################################################
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,3> >& ray,const TRIANGULATED_SURFACE<T>& surface, const T thickness_over_two)
{
    assert(surface.hierarchy);assert(surface.bounding_box);assert(surface.triangle_list);
    bool hit=false;T start_t,end_t;if(!INTERSECTION::Get_Intersection_Range(ray,*surface.bounding_box,start_t,end_t))return false;
    RANGE<VECTOR<T,3> > ray_box(ray.Point(start_t));ray_box.Enlarge_To_Include_Point(ray.Point(end_t));
    ARRAY<int> triangles_to_check;surface.hierarchy->Intersection_List(ray_box,triangles_to_check,4*thickness_over_two);
    for(int i=1;i<=triangles_to_check.m;i++){int k=triangles_to_check(i);if(INTERSECTION::Closest_Non_Intersecting_Point(ray,(*surface.triangle_list)(k),thickness_over_two)){ray.aggregate_id=k;hit=true;}}
    return hit;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const TRIANGULATED_SURFACE<float>&, const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,3> >&,const TRIANGULATED_SURFACE<float>&, const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,3> >&,const TRIANGULATED_SURFACE<double>&, const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,3> >&,const TRIANGULATED_SURFACE<double>&, const double);
#endif
};
};
