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
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Continuous_Collision_Detection/EDGE_EDGE_COLLISION.h>
using namespace PhysBAM;
using namespace CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS;
//#####################################################################
// Function Edge_Edge_Collision
//#####################################################################
template<class T,class T_ARRAY> bool PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Edge_Edge_Collision(const POINT_2D<T>& pt,const POINT_2D<T>& point,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V_edges,const T dt,const T collision_thickness,T& collision_time,
    VECTOR<T,2>& normal,VECTOR<T,2>& weights,T& relative_speed,bool allow_negative_weights,const T small_number,const bool exit_early)
{
    VECTOR<T,2> x2_minus_x1=point-pt,v2_minus_v1=V_edges(2)-V_edges(1);
    QUADRATIC<double> quadratic(v2_minus_v1.Magnitude_Squared(),2*VECTOR<T,2>::Dot_Product(x2_minus_x1,v2_minus_v1),x2_minus_x1.Magnitude_Squared()-sqr(collision_thickness));
    quadratic.Compute_Roots_In_Interval(0,dt);
    if(quadratic.roots==0)return false;
    else if(quadratic.roots==-1){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"VERY SINGULAR ON QUADRATIC SOLVE"<<std::endl;
#endif
        collision_time=0;}
    else if(quadratic.roots==1)collision_time=(T)quadratic.root1;
    else collision_time=min((T)quadratic.root1,(T)quadratic.root2);
    POINT_2D<T> new_point(pt+collision_time*V_edges(1));
    T distance;
    return new_point.Edge_Edge_Interaction(POINT_2D<T>(point+collision_time*V_edges(2)),V_edges,collision_thickness,distance,normal,weights,relative_speed,true,small_number,exit_early);
}
//#####################################################################
// Function Edge_Edge_Collision
//#####################################################################
template<class T> bool PhysBAM::CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::
Edge_Edge_Collision(const SEGMENT_3D<T>& seg_fault,const SEGMENT_3D<T>& segment,const VECTOR<T,3>& v1,const VECTOR<T,3>& v2,const VECTOR<T,3>& v3,const VECTOR<T,3>& v4,const T dt,
    const T collision_thickness,T& collision_time,VECTOR<T,3>& normal,VECTOR<T,2>& weights,T& relative_speed,const T small_number,const bool exit_early)
{
    // find cubic and compute the roots as possible collision times
    VECTOR<T,3> ABo=seg_fault.x2-seg_fault.x1,ABv=dt*(v2-v1),ACo=segment.x2-segment.x1,ACv=dt*(v4-v3);
    VECTOR<T,3> No=VECTOR<T,3>::Cross_Product(ABo,ACo),Nv=VECTOR<T,3>::Cross_Product(ABo,ACv)+VECTOR<T,3>::Cross_Product(ABv,ACo),Na=VECTOR<T,3>::Cross_Product(ABv,ACv);
    VECTOR<T,3> APo=segment.x1-seg_fault.x1,APv=dt*(v3-v1);
    
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
        T collision_time=dt*(T)iterative_solver.Bisection_Secant_Root(cubic,intervals(k).min_corner,intervals(k).max_corner);
        SEGMENT_3D<T> segment2(seg_fault.x1+collision_time*v1,seg_fault.x2+collision_time*v2);
        if(segment2.Edge_Edge_Interaction(SEGMENT_3D<T>(segment.x1+collision_time*v3,segment.x2+collision_time*v4),v1,v2,v3,v4,collision_thickness,distance,normal,weights,relative_speed,
                small_number,exit_early)) return true;}

    return false;
}
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Edge_Edge_Collision<float>(SEGMENT_3D<float> const&,SEGMENT_3D<float> const&,VECTOR<float,3> const&,VECTOR<float,3> const&,
    VECTOR<float,3> const&,VECTOR<float,3> const&,float,float,float&,VECTOR<float,3>&,VECTOR<float,2>&,float&,float,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Edge_Edge_Collision<float,ARRAY_VIEW<VECTOR<float,2> const,int> >(POINT_2D<float> const&,POINT_2D<float> const&,
    INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2> const,int>,VECTOR<int,2>&>,float,float,float&,VECTOR<float,2>&,VECTOR<float,2>&,float&,bool,float,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Edge_Edge_Collision<float,ARRAY_VIEW<VECTOR<float,2>,int> >(POINT_2D<float> const&,POINT_2D<float> const&,
    INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2>,int>,VECTOR<int,2>&>,float,float,float&,VECTOR<float,2>&,VECTOR<float,2>&,float&,bool,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Edge_Edge_Collision<double>(SEGMENT_3D<double> const&,SEGMENT_3D<double> const&,VECTOR<double,3> const&,VECTOR<double,3> const&,
    VECTOR<double,3> const&,VECTOR<double,3> const&,double,double,double&,VECTOR<double,3>&,VECTOR<double,2>&,double&,double,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Edge_Edge_Collision<double,ARRAY_VIEW<VECTOR<double,2> const,int> >(POINT_2D<double> const&,POINT_2D<double> const&,
    INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2> const,int>,VECTOR<int,2>&>,double,double,double&,VECTOR<double,2>&,VECTOR<double,2>&,double&,bool,double,bool);
template bool CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Edge_Edge_Collision<double,ARRAY_VIEW<VECTOR<double,2>,int> >(POINT_2D<double> const&,POINT_2D<double> const&,
    INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2>,int>,VECTOR<int,2>&>,double,double,double&,VECTOR<double,2>&,VECTOR<double,2>&,double&,bool,double,bool);
#endif
