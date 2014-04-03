//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_2D  
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Segment_Line_Intersection
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Segment_Line_Intersection(const VECTOR<T,2>& point_on_line,const VECTOR<T,2>& normal_of_line,T &interpolation_fraction) const
{
    T denominator=VECTOR<T,2>::Dot_Product(x2-x1,normal_of_line);
    if(!denominator){interpolation_fraction=FLT_MAX;return false;} // parallel
    interpolation_fraction=VECTOR<T,2>::Dot_Product(point_on_line-x1,normal_of_line)/denominator;
    return interpolation_fraction<=1 && interpolation_fraction>=0;
}
//#####################################################################
// Function Closest_Point_On_Segment
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_2D<T>::
Closest_Point_On_Segment(const VECTOR<T,2>& point) const
{                  
    VECTOR<T,2> v=x2-x1;
    T denominator=VECTOR<T,2>::Dot_Product(v,v);
    if(denominator == 0) return x1; // x1 and x2 are a single point
    else{
        T t=VECTOR<T,2>::Dot_Product(point-x1,v)/denominator;
        if(t <= 0) return x1;
        else if(t >= 1) return x2;
        else{v=x1+(x2-x1)*t;return v;}}
}
//#####################################################################
// Function Distance_From_Point_To_Segment
//#####################################################################
template<class T> T SEGMENT_2D<T>::
Distance_From_Point_To_Segment(const VECTOR<T,2>& point) const
{                  
    VECTOR<T,2> v=Closest_Point_On_Segment(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Closest_Point_On_Line
//#####################################################################
template<class T> VECTOR<T,2> SEGMENT_2D<T>::
Closest_Point_On_Line(const VECTOR<T,2>& point) const
{                  
    VECTOR<T,2> v=x2-x1;
    T denominator=VECTOR<T,2>::Dot_Product(v,v);
    if(denominator == 0) return x1; // x1 and x2 are a single point
    else{
        T t=VECTOR<T,2>::Dot_Product(point-x1,v)/denominator;
        v=x1+(x2-x1)*t;return v;}
}
//#####################################################################
// Function Distance_From_Point_To_Line
//#####################################################################
template<class T> T SEGMENT_2D<T>::
Distance_From_Point_To_Line(const VECTOR<T,2>& point) const
{                  
    VECTOR<T,2> v=Closest_Point_On_Line(point),d=v-point;
    return d.Magnitude();
}
//#####################################################################
// Function Shortest_Vector_Between_Segments
//#####################################################################
// vector points from input segment to this segment
// not accurate as the segments become parallel
template<class T> VECTOR<T,2> SEGMENT_2D<T>::
Shortest_Vector_Between_Segments(const SEGMENT_2D<T>& segment,T& a,T& b) const
{
    VECTOR<T,2> u=x2-x1,v=segment.x2-segment.x1,w=segment.x1-x1;
    T u_magnitude=u.Magnitude(),v_magnitude=v.Magnitude();
    T u_dot_u=sqr(u_magnitude),v_dot_v=sqr(v_magnitude),u_dot_v=VECTOR<T,2>::Dot_Product(u,v),
               u_dot_w=VECTOR<T,2>::Dot_Product(u,w),v_dot_w=VECTOR<T,2>::Dot_Product(v,w);
    T denominator=u_dot_u*v_dot_v-sqr(u_dot_v);
    T rhs1=v_dot_v*u_dot_w-u_dot_v*v_dot_w,rhs2=u_dot_v*u_dot_w-u_dot_u*v_dot_w;
    if(rhs1 <= 0) a=0;else if(denominator < rhs1) a=1;else a=rhs1/denominator;
    if(rhs2 <= 0) b=0;else if(denominator < rhs2) b=1;else b=rhs2/denominator;
    if(a > 0 && a < 1){
        if(b == 0){a=u_dot_w/u_dot_u;if(a < 0) a=0;else if(a > 1) a=1;}
        else if(b == 1){a=VECTOR<T,2>::Dot_Product(segment.x2-x1,u)/u_dot_u;if(a < 0) a=0;else if(a > 1) a=1;}}
    else if(b > 0 && b < 1){
        if(a == 0){b=-v_dot_w/v_dot_v;if(b < 0) b=0;else if(b > 1) b=1;}
        else if(a == 1){b=VECTOR<T,2>::Dot_Product(x2-segment.x1,v)/v_dot_v;if(b < 0) b=0;else if(b > 1) b=1;}}
    else{
        T a_distance=(a == 0 ? -rhs1:rhs1-denominator)*u_magnitude,b_distance=(b == 0 ? -rhs2:rhs2-denominator)*v_magnitude;
        if(a_distance > b_distance){
            if(a == 0) b=-v_dot_w/v_dot_v;else b=VECTOR<T,2>::Dot_Product(x2-segment.x1,v)/v_dot_v;if(b < 0) b=0;else if(b > 1) b=1;}
        else{if(b == 0) a=u_dot_w/u_dot_u;else a=VECTOR<T,2>::Dot_Product(segment.x2-x1,u)/u_dot_u;if(a < 0) a=0;else if(a > 1) a=1;}}
    return a*u-w-b*v;
}
//#####################################################################
// Function Segment_Segment_Interaction
//#####################################################################
template<class T> int SEGMENT_2D<T>::
Segment_Segment_Interaction(const SEGMENT_2D<T>& segment,const VECTOR<T,2>& v1,const VECTOR<T,2>& v2,const VECTOR<T,2>& v3,const VECTOR<T,2>& v4,const T interaction_distance,
                            T& distance,VECTOR<T,2>& normal,T& a,T& b,T& relative_speed,const T small_number) const
{
    normal=Shortest_Vector_Between_Segments(segment,a,b);
    distance=normal.Magnitude();if(distance > interaction_distance) return 0; // no interaction
    VECTOR<T,2> velocity1=(1-a)*v1+a*v2,velocity2=(1-b)*v3+b*v4;
    if(distance > small_number) normal/=distance;
    else{ // set normal based on relative velocity perpendicular to the two points
        VECTOR<T,2> relative_velocity=velocity1-velocity2;
        VECTOR<T,2> u=x2-x1;
        normal=relative_velocity-VECTOR<T,2>::Dot_Product(relative_velocity,u)/VECTOR<T,2>::Dot_Product(u,u)*u;
        T normal_magnitude=normal.Magnitude();
        if(normal_magnitude > small_number) normal/=normal_magnitude;
        else{ // relative velocity perpendicular to the segment is 0, pick any direction perpendicular to the segment
            if(abs(u.x) > abs(u.y) && abs(u.x)) normal=VECTOR<T,2>(0,1);
            else if(abs(u.y) > abs(u.x) && abs(u.y)) normal=VECTOR<T,2>(1,0);
            else normal=VECTOR<T,2>(0,1);
            normal=normal-VECTOR<T,2>::Dot_Product(normal,u)/VECTOR<T,2>::Dot_Product(u,u)*u;normal.Normalize();
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            LOG::cout << "                                            PICKING RANDOM NORMAL !!!!!!!!!!!!!!!!!!!!!!!" <<  std::endl;
#endif
    }}
    relative_speed=VECTOR<T,2>::Dot_Product(velocity1-velocity2,normal); // relative speed is in the normal direction
    return 1;
}
//#####################################################################
// Function Segment_Segment_Collision
//#####################################################################
// Needs to be fixed
#if 0
template<class T> int SEGMENT_2D<T>::
Segment_Segment_Collision(const SEGMENT_2D<T>& segment,const VECTOR<T,2>& v1,const VECTOR<T,2>& v2,const VECTOR<T,2>& v3,const VECTOR<T,2>& v4,const T dt,
                          const T collision_thickness,T& collision_time,VECTOR<T,2>& normal,T& a,T& b,T& relative_speed,const T small_number) const
{
    // find cubic and compute the roots as possible collision times
    VECTOR<T,2> ABo=x2-x1,ABv=dt*(v2-v1),ACo=segment.x2-segment.x1,ACv=dt*(v4-v3);
    VECTOR<T,3> No_3D=VECTOR<T,2>::Cross_Product(ABo,ACo),
                             Nv_3D=VECTOR<T,2>::Cross_Product(ABo,ACv)+VECTOR<T,2>::Cross_Product(ABv,ACo),
                             Na_3D=VECTOR<T,2>::Cross_Product(ABv,ACv);
    VECTOR<T,2> No(No_3D.x,No_3D.y),Nv(Nv_3D.x,Nv_3D.y), Na(Na_3D.x,Na_3D.y);
    VECTOR<T,2> APo=segment.x1-x1,APv=dt*(v3-v1);
    CUBIC<T> cubic(VECTOR<T,2>::Dot_Product(Na,APv),VECTOR<T,2>::Dot_Product(Nv,APv)+VECTOR<T,2>::Dot_Product(Na,APo),
                              VECTOR<T,2>::Dot_Product(No,APv)+VECTOR<T,2>::Dot_Product(Nv,APo),VECTOR<T,2>::Dot_Product(No,APo));
    cubic.Compute_Roots_In_Interval(0,1+1e-6);

    // check the collision times
    T distance;
    for(int roots=1;roots<=cubic.roots;roots++){
        if(roots == 1) collision_time=dt*cubic.root1;else if(roots == 2) collision_time=dt*cubic.root2;else collision_time=dt*cubic.root3;
        SEGMENT_2D segment2(x1+collision_time*v1,x2+collision_time*v2);
        if(segment2.Segment_Segment_Interaction(SEGMENT_2D(segment.x1+collision_time*v3,segment.x2+collision_time*v4),v1,v2,v3,v4,collision_thickness,distance,normal,a,b,relative_speed,
            small_number)) return 1;}
    return 0;
}
#endif
//#####################################################################
// Function Thickened_Box 
//#####################################################################
template<class T> ORIENTED_BOX<VECTOR<T,2> > SEGMENT_2D<T>::
Thickened_Oriented_Box(const T thickness_over_two) const 
{
    // make norm and tangent direction of thickness_over_two length
    VECTOR<T,2> segment_vector=x2-x1,segment_tangent=segment_vector.Normalized()*thickness_over_two,segment_normal(-segment_tangent.y,segment_tangent.x);
    // Form box point and intersect
    return ORIENTED_BOX<TV>(x1-segment_tangent-segment_normal,segment_vector+segment_tangent*2,segment_normal*2);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Inside(const VECTOR<T,2>& point,const T thickness_over_two) const 
{
    ORIENTED_BOX<TV> thickened_oriented_box=Thickened_Oriented_Box(thickness_over_two);
    return thickened_oriented_box.Lazy_Inside(point);
}
//#####################################################################
// Function Linear_Point_Inside_Segment
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Linear_Point_Inside_Segment(const TV& X,const T thickness_over_2) const
{
    VECTOR<T,2> weights=Barycentric_Coordinates(X);
    return weights.x>=-thickness_over_2 && weights.y>=-thickness_over_2;
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
// outputs unsigned distance
template<class T> bool SEGMENT_2D<T>::
Point_Face_Interaction(const VECTOR<T,2>& x,const T interaction_distance,const bool allow_negative_weights,T& distance) const
{      
    distance=VECTOR<T,2>::Dot_Product(x-x1,Normal());
    return abs(distance)<=interaction_distance && Linear_Point_Inside_Segment(x,allow_negative_weights?interaction_distance:0);
}
//#####################################################################
// Function Point_Face_Interaction_Data 
//#####################################################################
template<class T> void SEGMENT_2D<T>::
Point_Face_Interaction_Data(const VECTOR<T,2>& x,T& distance,VECTOR<T,2>& interaction_normal,VECTOR<T,2>& weights,const bool perform_attractions) const
{
    interaction_normal=Normal();weights=Barycentric_Coordinates(x);
    if(!perform_attractions && distance<0){distance*=-1;interaction_normal*=-1;} // distance > 0, interaction_normal points from the triangle to the point
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Point_Face_Interaction(const VECTOR<T,2>& x,const VECTOR<T,2>& v,const TV& v1,const TV& v2,const T interaction_distance,T& distance,
    VECTOR<T,2>& interaction_normal,VECTOR<T,2>& weights,T& relative_speed,const bool allow_negative_weights,const bool exit_early) const
{
    if(!Point_Face_Interaction(x,interaction_distance,allow_negative_weights,distance)) return false;
    if(!exit_early){
        Point_Face_Interaction_Data(x,distance,interaction_normal,weights,false);
        relative_speed=VECTOR<T,2>::Dot_Product(v-(weights.x*v1+weights.y*v2),interaction_normal);} // relative speed is in the normal direction
    return true;
}
//#####################################################################
// Function Robust_Point_Segment_Collision
//#####################################################################
template<class T> POINT_SIMPLEX_COLLISION_TYPE SEGMENT_2D<T>::
Robust_Point_Segment_Collision(const SEGMENT_2D<T>& initial_segment,const SEGMENT_2D<T>& final_segment,const VECTOR<T,2> &x,const VECTOR<T,2> &final_x,const T dt, const T collision_thickness,
    T& collision_time,TV& normal,T& collision_alpha,T& relative_speed)
{
    return CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Segment_Collision(initial_segment,final_segment,x,final_x,dt,collision_thickness,collision_time,normal,collision_alpha,
        relative_speed);
}
//#####################################################################
// Function Point_Face_Collision
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Point_Face_Collision(const TV& x,const TV& v,const TV& v1,const TV& v2,const T dt,const T collision_thickness,T& collision_time,TV& normal,VECTOR<T,2>& weights,T& relative_speed,
    const bool exit_early) const 
{
    return CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Face_Collision(*this,x,v,v1,v2,dt,collision_thickness,collision_time,normal,weights,relative_speed,exit_early);
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void SEGMENT_2D<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<SEGMENT_2D<T> >& clipped_simplices) const
{
    // cut with all sides of box
    clipped_simplices.Remove_All();
    clipped_simplices.Append(*this);
    for(int axis=1;axis<=TV::dimension;axis++){
        for(int i=clipped_simplices.m;i>=1;i--){
            SEGMENT_2D<T> triangle=clipped_simplices(i);clipped_simplices.Remove_Index_Lazy(i);
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(triangle,LINE_2D<T>(-TV::Axis_Vector(axis),box.min_corner),clipped_simplices);}
        for(int i=clipped_simplices.m;i>=1;i--){
            SEGMENT_2D<T> triangle=clipped_simplices(i);clipped_simplices.Remove_Index_Lazy(i);
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(triangle,LINE_2D<T>(TV::Axis_Vector(axis),box.max_corner),clipped_simplices);}}
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void SEGMENT_2D<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const SEGMENT_2D<T>& segment,const LINE_2D<T>& cutting_plane,ARRAY<SEGMENT_2D<T> >& negative_segments)
{
    VECTOR<T,2> phi_nodes;
    VECTOR<VECTOR<T,2>,2> X_nodes;
    X_nodes(1)=segment.x1;
    X_nodes(2)=segment.x2;
    for(int i=1;i<=2;i++){phi_nodes[i]=cutting_plane.Signed_Distance(X_nodes[i]);}
    int positive_count=0;
    for(int i=1;i<=2;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
        case 0: // in negative halfspace
            negative_segments.Append(segment);break;
        case 1:{
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            TV interface_location=LINEAR_INTERPOLATION<T,VECTOR<T,2> >::Linear(X_nodes[1],X_nodes[2],LEVELSET_UTILITIES<T>::Theta(phi_nodes[1],phi_nodes[2]));
            if(phi_nodes[1]>0) negative_segments.Append(SEGMENT_2D<T>(interface_location,segment.x2));
            else negative_segments.Append(SEGMENT_2D<T>(segment.x1,interface_location));
            break;}
        case 2: // in positive halfspace
            break;}
}
//#####################################################################
// Function Clip_To_Box_Helper
//#####################################################################
template<class T> bool
Clip_To_Box_Helper(T z1,T z2,T& a,T& b)
{
    if(z1>z2){z1=1-z1;z2=1-z2;}
    if(z1<0){
        if(z2<0) return false;
        a=max(a,z1/(z1-z2));}
    if(z2>1){
        if(z1>1) return false;
        b=min(b,(z1-1)/(z1-z2));}
    return true;
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> bool SEGMENT_2D<T>::
Clip_To_Box(const RANGE<TV>& box,T& a,T& b) const
{
    TV z1=(x1-box.min_corner)/box.Edge_Lengths();
    TV z2=(x2-box.min_corner)/box.Edge_Lengths();
    a=0;
    b=1;
    return Clip_To_Box_Helper(z1.x,z2.x,a,b) && Clip_To_Box_Helper(z1.y,z2.y,a,b) && a<b;
}
//#####################################################################
template class SEGMENT_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEGMENT_2D<double>;
#else
template bool SEGMENT_2D<double>::Point_Face_Interaction(VECTOR<double,2> const&,VECTOR<double,2> const&,VECTOR<double,2> const&,VECTOR<double,2> const&,double,double&,VECTOR<double,2>&,
    VECTOR<double,2>&,double&,bool,bool) const;
#endif
