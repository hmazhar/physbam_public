//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_ELLIPSOID_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const ELLIPSOID<T>& ellipsoid,const T thickness)
{ // TODO(jontg): everything else takes thickness_over_two; this should be made consistent...
    typedef VECTOR<T,3> TV;
    // Apply inverse transformation to ray
    TV ray_endpoint=ellipsoid.radii.Solve_Linear_System(ellipsoid.orientation.Inverse_Rotate(ray.endpoint-ellipsoid.center));
    TV ray_direction=ellipsoid.radii.Solve_Linear_System(ellipsoid.orientation.Inverse_Rotate(ray.direction)).Normalized();
    
    // sphere intersection routine
    T thickness_over_two=(T).5*thickness;
    T distance_squared=TV::Dot_Product(ray_endpoint,ray_endpoint);
    T outside_shell_squared=sqr((T)1+thickness_over_two);
    T c=distance_squared-outside_shell_squared;
    if(c > 0){ // outside   
        T b=TV::Dot_Product(ray_direction,ray_endpoint);if(b >= 0) return false; // no intersection - ray goes away from sphere 
        T d=sqr(b)-c;if(d < 0) return 0; // no intersection - ray misses sphere
        T t=-b-sqrt(d); // smaller of the two roots
        if(ray.semi_infinite || t < ray.t_max){ray.semi_infinite=false;ray.t_max=t;return true;}
        else return false;}
    else{ // inside or on the boundary
        T inside_shell_squared=sqr((T)1-thickness_over_two);
        T c=distance_squared-inside_shell_squared;
        if(c < 0){ // inside
            T b=TV::Dot_Product(ray_direction,ray_endpoint);
            T d=sqr(b)-c;
            T t=-b+sqrt(d); // larger of the two roots
            if(ray.semi_infinite || t < ray.t_max){ray.semi_infinite=false;ray.t_max=t;return true;}
            else return false;}
        else{ray.semi_infinite=false;ray.t_max=0;return true;}} // on the boundary
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const ELLIPSOID<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,3> >&,const ELLIPSOID<double>&,const double);
#endif
};
};
