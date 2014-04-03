//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_RING_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const RING<T>& ring,const T thickness)
{ // TODO(jontg): everything else takes thickness_over_two; this should be made consistent...
    typedef VECTOR<T,3> TV;
    bool intersection=false;T thickness_over_two=(T).5*thickness;
    
    // project ray into plane1 with x1 at the origin
    TV ray_endpoint=ray.endpoint-ring.plane1.x1;
    T vertical_endpoint=TV::Dot_Product(ray_endpoint,ring.plane1.normal);
    ray_endpoint-=vertical_endpoint*ring.plane1.normal;
    T vertical_direction=TV::Dot_Product(ray.direction,ring.plane1.normal);
    TV ray_direction=ray.direction-vertical_direction*ring.plane1.normal;

    // check intersections with the cylinder
    T distance_squared=TV::Dot_Product(ray_endpoint,ray_endpoint);
    T outer_shell_outside_squared=sqr(ring.outer_radius+thickness_over_two),outer_shell_inside_squared=sqr(ring.outer_radius-thickness_over_two);
    T inner_shell_outside_squared=sqr(ring.inner_radius-thickness_over_two),inner_shell_inside_squared=sqr(ring.inner_radius+thickness_over_two);
    if(distance_squared > outer_shell_outside_squared || distance_squared < inner_shell_outside_squared) {
        T b=TV::Dot_Product(ray_direction,ray_endpoint);
        T a=TV::Dot_Product(ray_direction,ray_direction); // note that ray_direction is not a unit vector
        // intersect with outer radius
        if(distance_squared > outer_shell_outside_squared){
            T c=distance_squared-outer_shell_outside_squared;
            if(b < 0){ // otherwise no intersection - ray goes away from cylinder
                T d=sqr(b)-a*c;
                if(d >= 0){ // otherwise no intersection - ray misses cylinder
                    T numerator;
                    numerator =-b-sqrt(d); // smaller of the two roots
                    if(a*vertical_endpoint+numerator*vertical_direction < a*thickness && 
                        a*vertical_endpoint+numerator*vertical_direction > -a*(ring.height+thickness)){ // intersects the finite cylinder
                        T t=numerator/a;
                        if(ray.semi_infinite || t < ray.t_max){intersection=true;ray.semi_infinite=false;ray.t_max=t;ray.aggregate_id=1;}}}}}
        // intersect with inner
        {T c=distance_squared-inner_shell_outside_squared;
        T d=sqr(b)-a*c;
        if(d >= 0){ // otherwise no intersection - ray misses cylinder
            T numerator=-b+sqrt(d); // larger of the two roots
            if(numerator>0 && a*vertical_endpoint+numerator*vertical_direction < a*thickness && a*vertical_endpoint+numerator*vertical_direction > -a*(ring.height+thickness)){ // intersects the finite cylinder
                T t=numerator/a;
                if(ray.semi_infinite || t < ray.t_max){intersection=true;ray.semi_infinite=false;ray.t_max=t;ray.aggregate_id=2;}}}}}
    else if(distance_squared < outer_shell_inside_squared && distance_squared > inner_shell_inside_squared) { // inside or on the boundary
        T b=TV::Dot_Product(ray_direction,ray_endpoint);
        T a=TV::Dot_Product(ray_direction,ray_direction); // note that ray_direction is not a unit vector
        // intersect with outer radius
        {T c=distance_squared-outer_shell_inside_squared;
        T d=sqr(b)-a*c;
        T numerator=-b+sqrt(d); // larger of the two roots
        if(a*vertical_endpoint+numerator*vertical_direction < a*thickness && a*vertical_endpoint+numerator*vertical_direction > -a*(ring.height+thickness)) { // intersects the finite cylinder
            T t=numerator/a;
            if(ray.semi_infinite || t < ray.t_max){intersection=true;ray.semi_infinite=false;ray.t_max=t;ray.aggregate_id=1;}}}
        // intersect with inner radius
        {T c=distance_squared-inner_shell_inside_squared;
         T d=sqr(b)-a*c;
         T numerator=-b-sqrt(d); // smaller of the two roots
         if(numerator>0 && a*vertical_endpoint+numerator*vertical_direction < a*thickness && a*vertical_endpoint+numerator*vertical_direction > -a*(ring.height+thickness)) { // intersects the finite cylinder
            T t=numerator/a;
            if(ray.semi_infinite || t < ray.t_max){intersection=true;ray.semi_infinite=false;ray.t_max=t;ray.aggregate_id=2;}}}}
    else{ // on the boundary
        if(vertical_endpoint < thickness_over_two && vertical_endpoint > -(ring.height+thickness_over_two)){ // intersects the finite cylinder
            intersection=true;ray.semi_infinite=false;ray.t_max=0;
            if(distance_squared >= outer_shell_inside_squared) ray.aggregate_id=1;
            else ray.aggregate_id=2;}}

    // check intersection with planes
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;
    if(INTERSECTION::Intersects(ray,ring.plane1,thickness)){
        TV point=ray.Point(ray.t_max);
        point-=TV::Dot_Product(point-ring.plane1.x1,ring.plane1.normal)*ring.plane1.normal;
        T distance=(point-ring.plane1.x1).Magnitude();
        if(distance <= ring.outer_radius+thickness && distance >= ring.inner_radius-thickness){intersection=true;ray.aggregate_id=3;}
        else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;}}
    else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;}

    save_semi_infinite=ray.semi_infinite;save_t_max=ray.t_max;save_aggregate=ray.aggregate_id;
    if(INTERSECTION::Intersects(ray,ring.plane2,thickness)){ // TODO(jontg) ... thickness?
        TV point=ray.Point(ray.t_max);
        point-=TV::Dot_Product(point-ring.plane1.x1,ring.plane1.normal)*ring.plane1.normal;
        T distance = (point-ring.plane1.x1).Magnitude();
        if(distance <= ring.outer_radius+thickness && distance >= ring.inner_radius-thickness){intersection=true;ray.aggregate_id=4;}
        else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;}}
    else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;}

    return intersection;

}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const RING<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,3> >&,const RING<double>&,const double);
#endif
};
};
