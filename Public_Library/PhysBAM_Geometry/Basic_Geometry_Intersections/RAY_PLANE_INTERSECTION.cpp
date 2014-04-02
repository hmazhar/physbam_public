//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane, const T thickness_over_two,const T distance,const T rate_of_approach)
{
    if(distance>-thickness_over_two && distance<thickness_over_two){ray.semi_infinite=false;ray.t_max=0;ray.intersection_location=RAY<VECTOR<T,3> >::START_POINT;return true;} // within the boundary
    if(rate_of_approach*distance<=0) return false; // no intersection
    if(!ray.semi_infinite){ // t_max is defined
        if(rate_of_approach>0 && distance-thickness_over_two<ray.t_max*rate_of_approach){
            ray.t_max=(distance-thickness_over_two)/rate_of_approach;ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
        else if(rate_of_approach<0 && distance+thickness_over_two>ray.t_max*rate_of_approach){
            ray.t_max=(distance+thickness_over_two)/rate_of_approach;ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
        else return false;} 
    else{ // semi-infinite ray with t_max=0 not defined
        if(rate_of_approach>0) ray.t_max=(distance-thickness_over_two)/rate_of_approach;
        else ray.t_max=(distance+thickness_over_two)/rate_of_approach;
        if(ray.t_max>FLT_MAX || ray.t_max<-FLT_MAX){ray.t_max=0;return false;} // guard against overflow
        ray.semi_infinite=false;ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane, const T thickness_over_two)
{
    T distance=VECTOR<T,3>::Dot_Product(plane.normal,ray.endpoint-plane.x1);
    if(distance>-thickness_over_two && distance<thickness_over_two){ray.semi_infinite=false;ray.t_max=0;ray.intersection_location=RAY<VECTOR<T,3> >::START_POINT;return true;} // within the boundary
    T rate_of_approach=-VECTOR<T,3>::Dot_Product(plane.normal,ray.direction);
    if(rate_of_approach*distance<=0) return false; // no intersection
    if(!ray.semi_infinite){ // t_max is defined
        if(rate_of_approach>0 && distance-thickness_over_two<ray.t_max*rate_of_approach){
            ray.t_max=(distance-thickness_over_two)/rate_of_approach;ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
        else if(rate_of_approach<0 && distance+thickness_over_two>ray.t_max*rate_of_approach){
            ray.t_max=(distance+thickness_over_two)/rate_of_approach;ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
        else return false;} 
    else{ // semi-infinite ray with t_max=0 not defined
        if(rate_of_approach>0) ray.t_max=(distance-thickness_over_two)/rate_of_approach;
        else ray.t_max=(distance+thickness_over_two)/rate_of_approach;
        if(ray.t_max>FLT_MAX || ray.t_max<-FLT_MAX){ray.t_max=0;return false;} // guard against overflow
        ray.semi_infinite=false;ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Lazy_Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane)
{
    T distance=VECTOR<T,3>::Dot_Product(plane.normal,ray.endpoint-plane.x1),rate_of_approach=-VECTOR<T,3>::Dot_Product(plane.normal,ray.direction);
    if(rate_of_approach*distance<=0) return false; // no intersection
    if(!ray.semi_infinite){ // t_max is defined
        if(rate_of_approach>0 && distance<ray.t_max*rate_of_approach){ray.t_max=distance/rate_of_approach;return true;}
        else if(rate_of_approach<0 && distance>ray.t_max*rate_of_approach){ray.t_max=distance/rate_of_approach;return true;}
        else return false;}
    else{ // semi-infinite ray with t_max=0 not defined
        ray.t_max=distance/rate_of_approach;
        if(ray.t_max>FLT_MAX || ray.t_max<-FLT_MAX){ray.t_max=0;return false;} // guard against overflow
        ray.semi_infinite=false;
        return true;}
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Rectangle_Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane,const PLANE<T>& bounding_plane1,const PLANE<T>& bounding_plane2,const PLANE<T>& bounding_plane3,const PLANE<T>& bounding_plane4,const T thickness_over_two)
{
    RAY<VECTOR<T,3> > ray_temp;ray.Save_Intersection_Information(ray_temp);
    if(INTERSECTION::Intersects(ray,plane,thickness_over_two)){
        VECTOR<T,3> point=ray.Point(ray.t_max);
        if(!bounding_plane1.Outside(point,thickness_over_two) && !bounding_plane2.Outside(point,thickness_over_two) && !bounding_plane3.Outside(point,thickness_over_two)
            && !bounding_plane4.Outside(point,thickness_over_two)) return true;
        else ray.Restore_Intersection_Information(ray_temp);}
    return false;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const PLANE<float>&,const float,const float,const float);
template bool Intersects(RAY<VECTOR<float,3> >&,const PLANE<float>&,const float);
template bool Lazy_Intersects(RAY<VECTOR<float,3> >&,const PLANE<float>&);
template bool Rectangle_Intersects(RAY<VECTOR<float,3> >&,const PLANE<float>&,const PLANE<float>&,const PLANE<float>&,const PLANE<float>&,const PLANE<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,3> >&,const PLANE<double>&,const double,const double,const double);
template bool Intersects(RAY<VECTOR<double,3> >&,const PLANE<double>&,const double);
template bool Lazy_Intersects(RAY<VECTOR<double,3> >&,const PLANE<double>&);
template bool Rectangle_Intersects(RAY<VECTOR<double,3> >&,const PLANE<double>&,const PLANE<double>&,const PLANE<double>&,const PLANE<double>&,const PLANE<double>&,const double);
#endif
};
};
