//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const RANGE<VECTOR<T,1> >& box,const T thickness_over_two)
{
    VECTOR<T,1> min_thickened=box.min_corner-thickness_over_two,max_thickened=box.max_corner+thickness_over_two;
    if(ray.direction.x>(T)0 && ray.endpoint.x<max_thickened.x){
        T distance_to_min=min_thickened.x-ray.endpoint.x;
        if(ray.semi_infinite || distance_to_min<=ray.t_max){
            ray.t_max=max((T)0,distance_to_min);ray.semi_infinite=false;return true;}
        else return false;}
    if(ray.direction.x<(T)0 && ray.endpoint.x>min_thickened.x){
        T distance_to_max=ray.endpoint.x-max_thickened.x;
        if(ray.semi_infinite || distance_to_max<=ray.t_max){
            ray.t_max=max((T)0,distance_to_max);ray.semi_infinite=false;return true;}
        else return false;}
    return false;
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,2> >& ray,const RANGE<VECTOR<T,2> >& box,const T thickness_over_two,const T segment_intersect_epsilon)
{
    VECTOR<T,2> min_thickened=box.min_corner-thickness_over_two,max_thickened=box.max_corner+thickness_over_two;
    bool left_test=  INTERSECTION::Intersection_Y_Segment(ray,min_thickened.x,min_thickened.y,max_thickened.y,segment_intersect_epsilon);
    bool right_test= INTERSECTION::Intersection_Y_Segment(ray,max_thickened.x,min_thickened.y,max_thickened.y,segment_intersect_epsilon);
    bool bottom_test=INTERSECTION::Intersection_X_Segment(ray,min_thickened.x,max_thickened.x,min_thickened.y,segment_intersect_epsilon);
    bool top_test=   INTERSECTION::Intersection_X_Segment(ray,min_thickened.x,max_thickened.x,max_thickened.y,segment_intersect_epsilon);
    return left_test || right_test || bottom_test || top_test;
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const RANGE<VECTOR<T,3> >& box,const T thickness_over_two)
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate_id=ray.aggregate_id,aggregate=0;
    VECTOR<T,3> min_thickened=box.min_corner-thickness_over_two,max_thickened=box.max_corner+thickness_over_two;
    VECTOR<T,3> corners[2]={box.min_corner,box.max_corner};

    for(int i=1;i<=3;++i){
        PLANE<T> plane(VECTOR<T,3>::Axis_Vector(i),VECTOR<T,3>());
        T rate_of_approach=-ray.direction(i);
        for(int j=0;j<=1;++j){
            plane.x1(i)=corners[j](i);
            if(Intersects(ray,plane,thickness_over_two,ray.endpoint(i)-corners[j](i),rate_of_approach)){
                VECTOR<T,2> point_no_i=ray.Point(ray.t_max).Remove_Index(i);
                if(point_no_i.All_Greater_Equal(min_thickened.Remove_Index(i)) && point_no_i.All_Less_Equal(max_thickened.Remove_Index(i))){
                    aggregate=save_aggregate_id=2*i+j-1;save_semi_infinite=false;save_t_max=ray.t_max;}
                else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate_id;}}}}

    if(aggregate){ray.aggregate_id=aggregate;return true;}
    return false;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool Lazy_Outside(RAY<VECTOR<T,3> >& ray,const RANGE<VECTOR<T,3> >& box)
{
    T t_min=0,local_t_max=ray.semi_infinite?FLT_MAX:ray.t_max;
    T new_t_min,inverse;
    #define Lazy_Outside_Dimension_Case(d,dmin,dmax,before,first_dimension)do{ \
        if(dmax before ray.endpoint.d)return false;inverse=1/ray.direction.d; \
        if(ray.endpoint.d before dmin){ \
            new_t_min=(ray.endpoint.d-dmin)*inverse; \
            t_min=first_dimension?new_t_min:max(t_min,new_t_min);} \
            local_t_max=min(local_t_max,(ray.endpoint.d-dmax)*inverse);if(local_t_max<t_min) return false;}while(0)
    if(ray.direction.x>0)Lazy_Outside_Dimension_Case(x,box.min_corner.x,box.max_corner.x,<,true);
    else if(ray.direction.x<0)Lazy_Outside_Dimension_Case(x,box.max_corner.x,box.min_corner.x,>,true);
    else if(ray.endpoint.x<box.min_corner.x || ray.endpoint.x>box.max_corner.x)return false;
    if(ray.direction.y>0)Lazy_Outside_Dimension_Case(y,box.min_corner.y,box.max_corner.y,<,false);
    else if(ray.direction.y<0)Lazy_Outside_Dimension_Case(y,box.max_corner.y,box.min_corner.y,>,false);
    else if(ray.endpoint.y<box.min_corner.y || ray.endpoint.y>box.max_corner.y)return false;
    if(ray.direction.z>0)Lazy_Outside_Dimension_Case(z,box.min_corner.z,box.max_corner.z,<,false);
    else if(ray.direction.z<0)Lazy_Outside_Dimension_Case(z,box.max_corner.z,box.min_corner.z,>,false);
    else if(ray.endpoint.z<box.min_corner.z || ray.endpoint.z>box.max_corner.z)return false;
    return true;
}
//#####################################################################
// Function Lazy_Intersects
//#####################################################################
template<class T> bool Lazy_Intersects(RAY<VECTOR<T,3> >& ray,const RANGE<VECTOR<T,3> >& box,const T thickness_over_two)
{
    // This comes from a paper "An efficient and Robust Ray-Box Intersection Algorithm" by Williams, Barrus, Morley, and Shirley
    // http://www.cs.utah.edu/~rmorley/pubs/box.pdf
    ray.Compute_Lazy_Box_Intersection_Acceleration_Data();
    VECTOR<T,3> extremes[2];
    extremes[0]=box.min_corner-thickness_over_two; // low extreme
    extremes[1]=box.max_corner+thickness_over_two; // high extreme
    T t_min=(extremes[ray.direction_is_negative.x].x-ray.endpoint.x)*ray.inverse_direction.x;
    T t_max=(extremes[1-ray.direction_is_negative.x].x-ray.endpoint.x)*ray.inverse_direction.x;
    T t_y_min=(extremes[ray.direction_is_negative.y].y-ray.endpoint.y)*ray.inverse_direction.y;
    T t_y_max=(extremes[1-ray.direction_is_negative.y].y-ray.endpoint.y)*ray.inverse_direction.y;
    if(t_min>t_y_max || t_y_min>t_max) return false;
    if(t_y_min>t_min)t_min=t_y_min;
    if(t_y_max<t_max)t_max=t_y_max;
    T t_z_min=(extremes[ray.direction_is_negative.z].z-ray.endpoint.z)*ray.inverse_direction.z;
    T t_z_max=(extremes[1-ray.direction_is_negative.z].z-ray.endpoint.z)*ray.inverse_direction.z;
    if(t_min>t_z_max || t_z_min>t_max) return false;
    if(t_z_min>t_min)t_min=t_z_min;
    if(t_z_max<t_max)t_max=t_z_max;
    //TODO:return aggregates properly
    if(t_max>0 && (ray.semi_infinite || t_min<ray.t_max)){ray.t_max=t_min;ray.semi_infinite=false;ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
    else return false;
}
//#####################################################################
// Function Get_Intersection_Range
//#####################################################################
template<class T> bool Get_Intersection_Range(const RAY<VECTOR<T,3> >& ray,const BOX<VECTOR<T,3> >& box,T& start_t,T& end_t)
{
    if(box.Lazy_Inside(ray.endpoint)) start_t=0;
    else{
        RAY<VECTOR<T,3> > ray_copy(ray);
        if(INTERSECTION::Lazy_Intersects(ray_copy,box)) start_t=ray_copy.t_max;
        else return false;}

    if(!ray.semi_infinite && box.Lazy_Inside(ray.Point(ray.t_max))) end_t=ray.t_max;
    else{
        RAY<VECTOR<T,3> > ray_copy(ray);
        if(ray_copy.semi_infinite){ray_copy.semi_infinite=false;ray_copy.t_max=start_t+2*box.Edge_Lengths().Max();}
        end_t=ray_copy.t_max; // save this
        ray_copy.endpoint=ray_copy.Point(ray_copy.t_max);ray_copy.direction*=-1; // flip ray
        if(INTERSECTION::Lazy_Intersects(ray_copy,box)) end_t-=ray_copy.t_max;
        else return false;}

    return true;
}
//#####################################################################
// Function Get_Intersection_Range
//#####################################################################
// template<class T> bool Get_Intersection_Range(const RAY<VECTOR<T,3> >& ray,const BOX<VECTOR<T,3> >& box,T& start_t,T& end_t)
template<class TV> bool Get_Intersection_Range(const RAY<TV>& ray,const RANGE<TV>& box,typename TV::SCALAR& start_t,typename TV::SCALAR& end_t)
{
    if(box.Lazy_Inside(ray.endpoint)) start_t=0;
    else{
        RAY<TV> ray_copy(ray);
        if(INTERSECTION::Intersects(ray_copy,box)) start_t=ray_copy.t_max;
        else return false;}

    if(!ray.semi_infinite && box.Lazy_Inside(ray.Point(ray.t_max))) end_t=ray.t_max;
    else{
        RAY<TV> ray_copy(ray);
        if(ray_copy.semi_infinite){ray_copy.semi_infinite=false;ray_copy.t_max=start_t+2*box.Edge_Lengths().Max();}
        end_t=ray_copy.t_max; // save this
        ray_copy.endpoint=ray_copy.Point(ray_copy.t_max);ray_copy.direction*=-1; // flip ray
        if(INTERSECTION::Intersects(ray_copy,box)) end_t-=ray_copy.t_max;
        else return false;}

    return true;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,1> >&,const RANGE<VECTOR<float,1> >&,const float);
template bool Intersects(RAY<VECTOR<float,2> >&,const RANGE<VECTOR<float,2> >&,const float,const float);
template bool Intersects(RAY<VECTOR<float,3> >&,const RANGE<VECTOR<float,3> >&,const float);
template bool Lazy_Outside(RAY<VECTOR<float,3> >&,const RANGE<VECTOR<float,3> >&);
template bool Lazy_Intersects(RAY<VECTOR<float,3> >&,const RANGE<VECTOR<float,3> >&,const float);
template bool Get_Intersection_Range(const RAY<VECTOR<float,3> >&,const BOX<VECTOR<float,3> >&,float&,float&);
template bool Get_Intersection_Range(const RAY<VECTOR<float,1> >&,const RANGE<VECTOR<float,1> >&,float&,float&);
template bool Get_Intersection_Range(const RAY<VECTOR<float,2> >&,const RANGE<VECTOR<float,2> >&,float&,float&);
template bool Get_Intersection_Range(const RAY<VECTOR<float,3> >&,const RANGE<VECTOR<float,3> >&,float&,float&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,1> >&,const RANGE<VECTOR<double,1> >&,const double);
template bool Intersects(RAY<VECTOR<double,2> >&,const RANGE<VECTOR<double,2> >&,const double,const double);
template bool Intersects(RAY<VECTOR<double,3> >&,const RANGE<VECTOR<double,3> >&,const double);
template bool Lazy_Outside(RAY<VECTOR<double,3> >&,const RANGE<VECTOR<double,3> >&);
template bool Lazy_Intersects(RAY<VECTOR<double,3> >&,const RANGE<VECTOR<double,3> >&,const double);
template bool Get_Intersection_Range(const RAY<VECTOR<double,3> >&,const BOX<VECTOR<double,3> >&,double&,double&);
template bool Get_Intersection_Range(const RAY<VECTOR<double,1> >&,const RANGE<VECTOR<double,1> >&,double&,double&);
template bool Get_Intersection_Range(const RAY<VECTOR<double,2> >&,const RANGE<VECTOR<double,2> >&,double&,double&);
template bool Get_Intersection_Range(const RAY<VECTOR<double,3> >&,const RANGE<VECTOR<double,3> >&,double&,double&);
#endif
};
};
