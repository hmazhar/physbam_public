//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_ORIENTED_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,2> >& ray,const SEGMENT_2D<T>& segment,const T thickness_over_two)
{
    VECTOR<T,2> from_start_to_start=segment.x1-ray.endpoint,segment_direction=segment.x2-segment.x1;T segment_length=segment_direction.Normalize();
    T cross_product=VECTOR<T,2>::Cross_Product(ray.direction,segment_direction).x,abs_cross_product=abs(cross_product);
    if(segment.Inside(ray.endpoint,thickness_over_two)){
        ray.t_max=0;ray.intersection_location=RAY<VECTOR<T,2> >::START_POINT;
        return true;}
    if((ray.semi_infinite && abs_cross_product>0) || (!ray.semi_infinite && ray.t_max*abs_cross_product>thickness_over_two)){
        T cross_recip=((T)1)/cross_product;
        T ray_t=cross_recip*VECTOR<T,2>::Cross_Product(from_start_to_start,segment_direction).x;
        if (ray_t<0||(ray_t>ray.t_max&&!ray.semi_infinite))return false;
        T segment_t=cross_recip*VECTOR<T,2>::Cross_Product(from_start_to_start,ray.direction).x;
        if (segment_t<-thickness_over_two||segment_t>segment_length+thickness_over_two)return false;
        ray.t_max=ray_t;ray.semi_infinite=false;return true;}
    return false;
}
//#####################################################################
// Function Fuzzy_Intersects
//#####################################################################
template<class T> bool Fuzzy_Intersects(RAY<VECTOR<T,2> >& ray,const SEGMENT_2D<T>& segment,const T thickness_over_two)
{
    // The order of these checks is important -- in particular we want to check 
    // intersection with lengthened segment before we try the endpoint in thickened box.
    ORIENTED_BOX<VECTOR<T,2> > thickened_oriented_box=segment.Thickened_Oriented_Box(thickness_over_two);
    if (thickened_oriented_box.Lazy_Inside(ray.endpoint)){
        ray.semi_infinite=false; ray.t_max=0;ray.intersection_location=RAY<VECTOR<T,2> >::START_POINT; return true;}
    else if (INTERSECTION::Intersects(ray,segment,thickness_over_two)){
        ray.intersection_location=RAY<VECTOR<T,2> >::INTERIOR_POINT; return true;}
    else if (!ray.semi_infinite && thickened_oriented_box.Lazy_Inside(ray.Point(ray.t_max))){
        ray.intersection_location=RAY<VECTOR<T,2> >::END_POINT; return true;}
    else
        return false;
}
//#####################################################################
// Function Closest_Non_Intersecting_Point
//#####################################################################
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,2> >& ray,const SEGMENT_2D<T>& segment,const T thickness_over_two)
{
    RAY<VECTOR<T,2> > ray_temp;ray.Save_Intersection_Information(ray_temp);
    if(!INTERSECTION::Fuzzy_Intersects(ray,segment,thickness_over_two)) return false;
    else if(ray.intersection_location==RAY<VECTOR<T,2> >::START_POINT) return true;
    else ray.Restore_Intersection_Information(ray_temp);

    // TODO: Save having to re-generate thickened oriented box (already generated in Fuzzy_Intersection)
    ORIENTED_BOX<VECTOR<T,2> > thickened_oriented_box=segment.Thickened_Oriented_Box(2*thickness_over_two);
    return INTERSECTION::Fuzzy_Intersects(ray,thickened_oriented_box,thickness_over_two);
}
//#####################################################################
// Function Intersection_X_Segment
//#####################################################################
// Optimized intersection for segment(x1,y),(x2,y), must have x1<x2
// Segment is lengthened at each end by thickness_over_two
template<class T> bool Intersection_X_Segment(RAY<VECTOR<T,2> >& ray,const T x1,const T x2,const T y,const T thickness_over_two)
{
    assert(x1<x2);
    VECTOR<T,2> from_start_to_start(x1-ray.endpoint.x,y-ray.endpoint.y);T length=x2-x1;
    T cross_product=-ray.direction.y,abs_cross_product=abs(cross_product);
    if((ray.semi_infinite && abs_cross_product>0) || (!ray.semi_infinite && ray.t_max*abs_cross_product>thickness_over_two)){
        T cross_recip=((T)1)/cross_product;
        T ray_t=-cross_recip*from_start_to_start.y;
        if (ray_t<0||(ray_t>ray.t_max&&!ray.semi_infinite))return false;
        T t=cross_recip*VECTOR<T,2>::Cross_Product(from_start_to_start,ray.direction).x;
        if (t<-thickness_over_two||t>length+thickness_over_two)return false;
        ray.t_max=ray_t;ray.semi_infinite=false;return true;}
    return false;
}
//#####################################################################
// Function Intersection_Y_Segment
//#####################################################################
// Optimized intersection for segment (x,y1),(x,y2), must have y1<y2
// Segment is lengthened at each end by thickness_over_two
template<class T> bool Intersection_Y_Segment(RAY<VECTOR<T,2> >& ray,const T x,const T y1,const T y2,const T thickness_over_two)
{
    assert(y1<y2);
    VECTOR<T,2> from_start_to_start(x-ray.endpoint.x,y1-ray.endpoint.y);T segment_length=y2-y1;
    T cross_product=ray.direction.x,abs_cross_product=abs(cross_product);
    if((ray.semi_infinite && abs_cross_product>0) || (!ray.semi_infinite && ray.t_max*abs_cross_product>thickness_over_two)){
        T cross_recip=((T)1)/cross_product;
        T ray_t=cross_recip*from_start_to_start.x;
        if (ray_t<0||(ray_t>ray.t_max&&!ray.semi_infinite))return false;
        T segment_t=cross_recip*VECTOR<T,2>::Cross_Product(from_start_to_start,ray.direction).x;
        if (segment_t<-thickness_over_two||segment_t>segment_length+thickness_over_two)return false;
        ray.t_max=ray_t;ray.semi_infinite=false;return true;}
    return false;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,2> >&,const SEGMENT_2D<float>&,const float);
template bool Fuzzy_Intersects(RAY<VECTOR<float,2> >&,const SEGMENT_2D<float>&,const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,2> >&,const SEGMENT_2D<float>&,const float);
template bool Intersection_X_Segment(RAY<VECTOR<float,2> >&,const float,const float,const float,const float);
template bool Intersection_Y_Segment(RAY<VECTOR<float,2> >&,const float,const float,const float,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,2> >&,const SEGMENT_2D<double>&,const double);
template bool Fuzzy_Intersects(RAY<VECTOR<double,2> >&,const SEGMENT_2D<double>&,const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,2> >&,const SEGMENT_2D<double>&,const double);
template bool Intersection_X_Segment(RAY<VECTOR<double,2> >&,const double,const double,const double,const double);
template bool Intersection_Y_Segment(RAY<VECTOR<double,2> >&,const double,const double,const double,const double);
#endif
};
};
