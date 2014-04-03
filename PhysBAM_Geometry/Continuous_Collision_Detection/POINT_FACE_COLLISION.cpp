//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.h>
using namespace PhysBAM;
using namespace CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS;
//#####################################################################
// Function Robust_Point_Triangle_Collision
//#####################################################################
template<class T> POINT_SIMPLEX_COLLISION_TYPE PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Robust_Point_Triangle_Collision(const TRIANGLE_3D<T>& initial_triangle,const TRIANGLE_3D<T>& final_triangle,const VECTOR<T,3>& x,const VECTOR<T,3>& final_x,const T dt,
    const T collision_thickness,T& collision_time,VECTOR<T,3>& normal,VECTOR<T,3>& weights,T& relative_speed)
{
    if(final_triangle.Point_Inside_Triangle(final_x,collision_thickness)){
        collision_time=dt;weights=final_triangle.Barycentric_Coordinates(final_x);normal=final_triangle.normal;return POINT_SIMPLEX_COLLISION_ENDS_INSIDE;}
    else if(initial_triangle.Point_Inside_Triangle(x,collision_thickness)){
        collision_time=0;weights=initial_triangle.Barycentric_Coordinates(x);normal=initial_triangle.normal;return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
    else{
        VECTOR<T,3> v1=(final_triangle.x1-initial_triangle.x1)/dt,v2=(final_triangle.x2-initial_triangle.x2)/dt,v3=(final_triangle.x3-initial_triangle.x3)/dt,v=(final_x-x)/dt;
        if(initial_triangle.Point_Face_Collision(x,v,v1,v2,v3,dt,collision_thickness,collision_time,normal,weights,relative_speed,false)) return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;
        else return POINT_SIMPLEX_NO_COLLISION;}
}
//#####################################################################
// Function Point_Triangle_Collision
//#####################################################################
template<class T> bool PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Point_Face_Collision(const TRIANGLE_3D<T>& tri,const VECTOR<T,3>& x,const VECTOR<T,3>& v,const VECTOR<T,3>& v1,const VECTOR<T,3>& v2,const VECTOR<T,3>& v3,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,3>& normal,VECTOR<T,3>& weights,T& relative_speed,const bool exit_early)
{
    // find cubic and compute the roots as possible collision times
    VECTOR<T,3> ABo=tri.x2-tri.x1,ABv=dt*(v2-v1),ACo=tri.x3-tri.x1,ACv=dt*(v3-v1);
    VECTOR<T,3> No=VECTOR<T,3>::Cross_Product(ABo,ACo),Nv=VECTOR<T,3>::Cross_Product(ABo,ACv)+VECTOR<T,3>::Cross_Product(ABv,ACo),Na=VECTOR<T,3>::Cross_Product(ABv,ACv);
    VECTOR<T,3> APo=x-tri.x1,APv=dt*(v-v1);

    CUBIC<double> cubic((double)VECTOR<T,3>::Dot_Product(Na,APv),(double)VECTOR<T,3>::Dot_Product(Nv,APv)+VECTOR<T,3>::Dot_Product(Na,APo),
                                       (double)VECTOR<T,3>::Dot_Product(No,APv)+VECTOR<T,3>::Dot_Product(Nv,APo),(double)VECTOR<T,3>::Dot_Product(No,APo));
    double xmin=0,xmax=1.000001;
    int num_intervals=0;VECTOR<INTERVAL<double>,3> intervals;
    cubic.Compute_Intervals(xmin,xmax,num_intervals,intervals(1),intervals(2),intervals(3));
    if(!num_intervals) return false;
  
    // find and check roots
    T distance;
    ITERATIVE_SOLVER<double> iterative_solver;iterative_solver.tolerance=1e-14;
    for(int k=1;k<=num_intervals;k++){
        collision_time=dt*(T)iterative_solver.Bisection_Secant_Root(cubic,intervals(k).min_corner,intervals(k).max_corner);
        TRIANGLE_3D<T> triangle(tri.x1+collision_time*v1,tri.x2+collision_time*v2,tri.x3+collision_time*v3);
        if(triangle.Point_Face_Interaction(x+collision_time*v,v,v1,v2,v3,collision_thickness,distance,normal,weights,relative_speed,true,exit_early)) return true;}

    return false;
}
//#####################################################################
// Function Robust_Point_Segment_Collision
//#####################################################################
template<class T> POINT_SIMPLEX_COLLISION_TYPE PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Robust_Point_Segment_Collision(const SEGMENT_2D<T>& initial_segment,const SEGMENT_2D<T>& final_segment,const VECTOR<T,2> &x,const VECTOR<T,2> &final_x,const T dt, const T collision_thickness,
    T& collision_time,VECTOR<T,2>& normal,T& collision_alpha,T& relative_speed)
{
    if(final_segment.Thickened_Oriented_Box(collision_thickness).Lazy_Inside(final_x)){
        collision_time=dt;collision_alpha=final_segment.Barycentric_Coordinates(final_x).y;return POINT_SIMPLEX_COLLISION_ENDS_INSIDE;}
    else if(initial_segment.Thickened_Oriented_Box(collision_thickness).Lazy_Inside(x)){
        collision_time=0;collision_alpha=initial_segment.Barycentric_Coordinates(x).y;return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
    else{
        VECTOR<T,2> v1=(final_segment.x1-initial_segment.x1)/dt,v2=(final_segment.x2-initial_segment.x2)/dt,v=(final_x-x)/dt;
        VECTOR<T,2> weights;
        if(initial_segment.Point_Face_Collision(x,v,v1,v2,dt,collision_thickness,collision_time,normal,weights,relative_speed,false)){
            collision_alpha=weights.y;return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
        else return POINT_SIMPLEX_NO_COLLISION;}
}
//#####################################################################
// Function Point_Face_Collision
//#####################################################################
template<class T> bool PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Point_Face_Collision(const SEGMENT_2D<T>& seg_fault,const VECTOR<T,2>& x,const VECTOR<T,2>& v,const VECTOR<T,2>& v1,const VECTOR<T,2>& v2,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,2>& normal,VECTOR<T,2>& weights,T& relative_speed,const bool exit_early)
{
    VECTOR<double,2> v_minus_v1(v-v1),v2_minus_v1(v2-v1),x_minus_x1(x-seg_fault.x1),x2_minus_x1(seg_fault.x2-seg_fault.x1);
    QUADRATIC<double> quadratic(VECTOR<double,2>::Cross_Product(v_minus_v1,v2_minus_v1).x,
        VECTOR<double,2>::Cross_Product(x_minus_x1,v2_minus_v1).x+VECTOR<double,2>::Cross_Product(v_minus_v1,x2_minus_x1).x,
        VECTOR<double,2>::Cross_Product(x_minus_x1,x2_minus_x1).x);

    double collision_time_temp(0),relative_speed_temp(0),distance(0);VECTOR<double,2> normal_temp,weights_temp;
    if(abs(1e3*quadratic.a)<abs(quadratic.b)){
        collision_time_temp=-quadratic.c/quadratic.b;
        if(collision_time_temp<0 || collision_time_temp>dt) return false;}
    else{
        quadratic.Compute_Roots_In_Interval(0,dt);
        if(quadratic.roots==0)return false;
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        else if(quadratic.roots==-1){LOG::cout<<"VERY SINGULAR ON QUADRATIC SOLVE"<<std::endl;collision_time_temp=0;}
#endif
        else if(quadratic.roots==1)collision_time_temp=quadratic.root1;
        else collision_time_temp=quadratic.root1;}
    SEGMENT_2D<double> segment((VECTOR<double,2>)seg_fault.x1+collision_time_temp*(VECTOR<double,2>)v1,(VECTOR<double,2>)seg_fault.x2+collision_time_temp*(VECTOR<double,2>)v2);
    bool interaction=segment.Point_Face_Interaction((VECTOR<double,2>)x+collision_time_temp*(VECTOR<double,2>)v,(VECTOR<double,2>)v,(VECTOR<double,2>)v1,(VECTOR<double,2>)v2,
        (double)collision_thickness,distance,normal_temp,weights_temp,relative_speed_temp,true,exit_early);
    collision_time=(T)collision_time_temp;relative_speed=(T)relative_speed_temp;normal=(VECTOR<T,2>)normal_temp;weights=(VECTOR<T,2>)weights_temp;
    return interaction;
}
//#####################################################################
// Function Robust_Point_Point_Collision
//#####################################################################
template<class T> POINT_SIMPLEX_COLLISION_TYPE PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Robust_Point_Point_Collision(const POINT_SIMPLEX_1D<T>& initial_simplex,const POINT_SIMPLEX_1D<T>& final_simplex,const VECTOR<T,1>& x,const VECTOR<T,1>& final_x,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,1>& normal,ONE& weights,T& relative_speed)
{
    if(final_simplex.Bounding_Box().Thickened(collision_thickness).Lazy_Inside(final_x)){
        collision_time=dt;return POINT_SIMPLEX_COLLISION_ENDS_INSIDE;}
    else if(initial_simplex.Bounding_Box().Thickened(collision_thickness).Lazy_Inside(x)){
        collision_time=0;return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
    else{
        VECTOR<T,1> v1=(final_simplex.x1-initial_simplex.x1)/dt,v=(final_x-x)/dt;
        if(Point_Point_Collision(initial_simplex,x,v,v1,dt,collision_thickness,collision_time,normal,weights,relative_speed,false))
            return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;
        else return POINT_SIMPLEX_NO_COLLISION;}
}
//#####################################################################
// Function Point_Point_Collision
//#####################################################################
template<class T> bool PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Point_Point_Collision(const POINT_SIMPLEX_1D<T>& initial_simplex,const VECTOR<T,1>& x,const VECTOR<T,1>& v,const VECTOR<T,1>& v1,const T dt,const T collision_thickness,
    T& collision_time,VECTOR<T,1>& normal,ONE& weights,T& relative_speed,const bool exit_early)
{
    T distance=x.x-initial_simplex.x1.x;
    relative_speed=v.x-v1.x;

    if(distance*relative_speed >= 0 || abs(distance) > abs(dt*relative_speed)) return false;
    collision_time=distance/relative_speed;
    if(distance > 0) normal=VECTOR<T,1>(1); else normal=VECTOR<T,1>(-1);

    return true;
}
//#####################################################################
template POINT_SIMPLEX_COLLISION_TYPE CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Point_Collision<float>(POINT_SIMPLEX_1D<float> const&,POINT_SIMPLEX_1D<float> const&,
    VECTOR<float,1> const&,VECTOR<float,1> const&,float,float,float&,VECTOR<float,1>&,ONE&,float&);
template POINT_SIMPLEX_COLLISION_TYPE CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Segment_Collision<float>(SEGMENT_2D<float> const&,SEGMENT_2D<float> const&,
    VECTOR<float,2> const&,VECTOR<float,2> const&,float,float,float&,VECTOR<float,2>&,float&,float&);
template POINT_SIMPLEX_COLLISION_TYPE CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Triangle_Collision<float>(TRIANGLE_3D<float> const&,TRIANGLE_3D<float> const&,
    VECTOR<float,3> const&,VECTOR<float,3> const&,float,float,float&,VECTOR<float,3>&,VECTOR<float,3>&,float&);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Point_Collision(POINT_SIMPLEX_1D<float> const&,VECTOR<float,1> const&,VECTOR<float,1> const&,VECTOR<float,1> const&,
    float,float,float&,VECTOR<float,1>&,ONE&,float&,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Face_Collision<float>(SEGMENT_2D<float> const&,VECTOR<float,2> const&,VECTOR<float,2> const&,VECTOR<float,2> const&,
    VECTOR<float,2> const&,float,float,float&,VECTOR<float,2>&,VECTOR<float,2>&,float&,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Face_Collision<float>(TRIANGLE_3D<float> const&,VECTOR<float,3> const&,VECTOR<float,3> const&,VECTOR<float,3> const&,
    VECTOR<float,3> const&,VECTOR<float,3> const&,float,float,float&,VECTOR<float,3>&,VECTOR<float,3>&,float&,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template POINT_SIMPLEX_COLLISION_TYPE CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Point_Collision<double>(POINT_SIMPLEX_1D<double> const&,POINT_SIMPLEX_1D<double> const&,
    VECTOR<double,1> const&,VECTOR<double,1> const&,double,double,double&,VECTOR<double,1>&,ONE&,double&);
template POINT_SIMPLEX_COLLISION_TYPE CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Segment_Collision<double>(SEGMENT_2D<double> const&,SEGMENT_2D<double> const&,
    VECTOR<double,2> const&,VECTOR<double,2> const&,double,double,double&,VECTOR<double,2>&,double&,double&);
template POINT_SIMPLEX_COLLISION_TYPE CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Triangle_Collision<double>(TRIANGLE_3D<double> const&,TRIANGLE_3D<double> const&,
    VECTOR<double,3> const&,VECTOR<double,3> const&,double,double,double&,VECTOR<double,3>&,VECTOR<double,3>&,double&);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Point_Collision(POINT_SIMPLEX_1D<double> const&,VECTOR<double,1> const&,VECTOR<double,1> const&,VECTOR<double,1> const&,
    double,double,double&,VECTOR<double,1>&,ONE&,double&,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Face_Collision<double>(SEGMENT_2D<double> const&,VECTOR<double,2> const&,VECTOR<double,2> const&,VECTOR<double,2> const&,
    VECTOR<double,2> const&,double,double,double&,VECTOR<double,2>&,VECTOR<double,2>&,double&,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Face_Collision<double>(TRIANGLE_3D<double> const&,VECTOR<double,3> const&,VECTOR<double,3> const&,VECTOR<double,3> const&,
    VECTOR<double,3> const&,VECTOR<double,3> const&,double,double,double&,VECTOR<double,3>&,VECTOR<double,3>&,double&,bool);
#endif
