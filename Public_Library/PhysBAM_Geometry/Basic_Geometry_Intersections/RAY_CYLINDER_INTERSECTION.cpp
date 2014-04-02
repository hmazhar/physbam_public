//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_CYLINDER_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const CYLINDER<T>& cylinder,const T thickness)
{
    typedef VECTOR<T,3> TV;
    bool intersection=false;
    T thickness_over_two=(T).5*thickness;
    
    // project ray into cylinder.plane1 with x1 at the origin
    TV ray_endpoint=ray.endpoint-cylinder.plane1.x1;
    T vertical_endpoint=TV::Dot_Product(ray_endpoint,cylinder.plane1.normal);
    ray_endpoint-=vertical_endpoint*cylinder.plane1.normal;
    T vertical_direction=TV::Dot_Product(ray.direction,cylinder.plane1.normal);
    TV ray_direction=ray.direction-vertical_direction*cylinder.plane1.normal;

    // check intersections with the cylinder
    T distance_squared=TV::Dot_Product(ray_endpoint,ray_endpoint);
    T outside_shell_squared=sqr(cylinder.radius+thickness_over_two);
    T c=distance_squared-outside_shell_squared;
    if(c > 0){ // outside   
        T b=TV::Dot_Product(ray_direction,ray_endpoint);
        if(b < 0){ // otherwise no intersection - ray goes away from cylinder
            T a=TV::Dot_Product(ray_direction,ray_direction); // note that ray_direction is not a unit vector
            T d=sqr(b)-a*c;
            if(d >= 0){ // otherwise no intersection - ray misses cylinder
                T numerator=-b-sqrt(d); // smaller of the two roots
                if(a*vertical_endpoint+numerator*vertical_direction < a*thickness && 
                    a*vertical_endpoint+numerator*vertical_direction > -a*(cylinder.height+thickness)){ // intersects the finite cylinder
                    T t=numerator/a;
                    if(ray.semi_infinite || t < ray.t_max){intersection=true;ray.semi_infinite=false;ray.t_max=t;ray.aggregate_id=1;}}}}}
    else{ // inside or on the boundary
        T inside_shell_squared=sqr(cylinder.radius-thickness_over_two);
        T c=distance_squared-inside_shell_squared;
        if(c < 0){ // inside
            T b=TV::Dot_Product(ray_direction,ray_endpoint);
            T a=TV::Dot_Product(ray_direction,ray_direction); // note that ray_direction is not a unit vector
            T d=sqr(b)-a*c;
            T numerator=-b+sqrt(d); // larger of the two roots
            if(a*vertical_endpoint+numerator*vertical_direction < a*thickness && 
                a*vertical_endpoint+numerator*vertical_direction > -a*(cylinder.height+thickness)){ // intersects the finite cylinder
                T t=numerator/a;
                if(ray.semi_infinite || t < ray.t_max){intersection=true;ray.semi_infinite=false;ray.t_max=t;ray.aggregate_id=1;}}}
        else{ // on the boundary 
            if(vertical_endpoint < thickness_over_two && vertical_endpoint > -(cylinder.height+thickness_over_two)){ // intersects the finite cylinder
                intersection=true;ray.semi_infinite=false;ray.t_max=0;ray.aggregate_id=1;}}}

    // check intersection with planes
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;
    if(INTERSECTION::Intersects(ray,cylinder.plane1,thickness)){
        TV point=ray.Point(ray.t_max);
        point-=TV::Dot_Product(point-cylinder.plane1.x1,cylinder.plane1.normal)*cylinder.plane1.normal;
        if((point-cylinder.plane1.x1).Magnitude() <= cylinder.radius+thickness){intersection=true;ray.aggregate_id=2;}
        else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;}}
    if(INTERSECTION::Intersects(ray,cylinder.plane2,thickness)){
        TV point=ray.Point(ray.t_max);
        point-=TV::Dot_Product(point-cylinder.plane1.x1,cylinder.plane1.normal)*cylinder.plane1.normal;
        if((point-cylinder.plane1.x1).Magnitude() <= cylinder.radius+thickness){intersection=true;ray.aggregate_id=3;}
        else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;}}

    return intersection;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const CYLINDER<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,3> >&,const CYLINDER<double>&,const double);
#endif
};
};
