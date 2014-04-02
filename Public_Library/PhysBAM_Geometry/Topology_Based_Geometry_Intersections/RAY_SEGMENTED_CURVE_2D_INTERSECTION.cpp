//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,2> >& ray,const SEGMENTED_CURVE_2D<T>& curve, const T thickness_over_two)
{
    if(curve.bounding_box){RAY<VECTOR<T,2> > tmp_ray=ray;
        if(curve.bounding_box->Outside(ray.endpoint,3*thickness_over_two) && !INTERSECTION::Intersects(tmp_ray,*curve.bounding_box,2*thickness_over_two,thickness_over_two)) return false;}
    bool hit=false;
    for(int i=1;i<=curve.mesh.elements.m;i++){
        SEGMENT_2D<T> segment=(curve.segment_list)?(*curve.segment_list)(i):curve.Get_Element(i);
        if(INTERSECTION::Fuzzy_Intersects(ray,segment,thickness_over_two)){ray.aggregate_id=i;hit=true;}}
    return hit;
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,2> >& ray,const SEGMENTED_CURVE_2D<T>& curve,const T thickness_over_two)
{
    if(curve.bounding_box){RAY<VECTOR<T,2> > tmp_ray=ray;
        if(curve.bounding_box->Outside(ray.endpoint,5*thickness_over_two) && !INTERSECTION::Intersects(tmp_ray,*curve.bounding_box,4*thickness_over_two,thickness_over_two)) return false;}
    bool hit=false;
    for(int i=1;i<=curve.mesh.elements.m;i++){
        SEGMENT_2D<T> segment=(curve.segment_list)?(*curve.segment_list)(i):curve.Get_Element(i);
        if(INTERSECTION::Closest_Non_Intersecting_Point(ray,segment,thickness_over_two)){ray.aggregate_id=i;hit=true;}}
    return hit;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,2> >&,const SEGMENTED_CURVE_2D<float>&, const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,2> >&,const SEGMENTED_CURVE_2D<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,2> >&,const SEGMENTED_CURVE_2D<double>&, const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,2> >&,const SEGMENTED_CURVE_2D<double>&,const double);
#endif
};
};
