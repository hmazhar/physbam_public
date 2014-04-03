//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SPHERE_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T,class TV> bool Intersects_helper(RAY<TV>& ray,const SPHERE<TV>& sphere,const T thickness)
{
    // translate ray so that the sphere can be considered at the origin
    TV ray_endpoint=ray.endpoint-sphere.center; 

    T thickness_over_two=(T).5*thickness;
    T distance_squared=ray_endpoint.Magnitude_Squared();
    T outside_shell_squared=sqr(sphere.radius+thickness_over_two);
    T c=distance_squared-outside_shell_squared;
    if(c > 0){ // outside   
        T b=TV::Dot_Product(ray.direction,ray_endpoint);if(b >= 0) return 0; // no intersection - ray goes away from sphere
        T d=sqr(b)-c;if(d < 0) return 0; // no intersection - ray misses sphere
        T t=-b-sqrt(d); // smaller of the two roots
        if(ray.semi_infinite || t < ray.t_max){ray.semi_infinite=false;ray.t_max=t;return 1;}
        else return 0;}
    else{ // inside or on the boundary
        T inside_shell_squared=sqr(sphere.radius-thickness_over_two);
        T c=distance_squared-inside_shell_squared;
        if(c < 0){ // inside
            T b=TV::Dot_Product(ray.direction,ray_endpoint);
            T d=sqr(b)-c;
            T t=-b+sqrt(d); // larger of the two roots
            if(ray.semi_infinite || t < ray.t_max){ray.semi_infinite=false;ray.t_max=t;return 1;}
            else return 0;}
        else{ray.semi_infinite=false;ray.t_max=0;return 1;}} // on the boundary
}
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const SPHERE<VECTOR<T,1> >& sphere,const T thickness){
    return Intersects_helper(ray,sphere,thickness);
}
template<class T> bool Intersects(RAY<VECTOR<T,2> >& ray,const SPHERE<VECTOR<T,2> >& sphere,const T thickness){
    return Intersects_helper(ray,sphere,thickness);
}
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const SPHERE<VECTOR<T,3> >& sphere,const T thickness){
    return Intersects_helper(ray,sphere,thickness);
}
//#####################################################################

template bool Intersects(RAY<VECTOR<float,1> >&,const SPHERE<VECTOR<float,1> >&,const float);
template bool Intersects(RAY<VECTOR<float,2> >&,const SPHERE<VECTOR<float,2> >&,const float);
template bool Intersects(RAY<VECTOR<float,3> >&,const SPHERE<VECTOR<float,3> >&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,1> >&,const SPHERE<VECTOR<double,1> >&,const double);
template bool Intersects(RAY<VECTOR<double,2> >&,const SPHERE<VECTOR<double,2> >&,const double);
template bool Intersects(RAY<VECTOR<double,3> >&,const SPHERE<VECTOR<double,3> >&,const double);
#endif
}
}
