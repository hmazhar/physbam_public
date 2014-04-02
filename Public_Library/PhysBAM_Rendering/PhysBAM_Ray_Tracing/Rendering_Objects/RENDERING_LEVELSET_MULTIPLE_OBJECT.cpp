//#####################################################################
// Copyright 2005, Jiayi Chong, Jeong-Mo Hong, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_ON_A_RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_LEVELSET_MULTIPLE_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Function Iterative_Solver_Tolerance
//#####################################################################
namespace{
template<class T> T Iterative_Solver_Tolerance(){STATIC_ASSERT((T)false);}
template<> float Iterative_Solver_Tolerance(){return (float).01;}
template<> double Iterative_Solver_Tolerance(){return .001;}
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T_LEVELSET_MULTIPLE> bool RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Intersection(RAY<VECTOR<T,3> >& ray,int &region_start,int &region_end,const T thickness) const
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;

    T t_start,t_end; 
    RANGE<VECTOR<T,3> > box=levelset_multiple.grid.domain;
    bool exit_intersection=false;int exit_aggregate=0;T exit_t_max=0;
    int intersect_box=INTERSECTION::Intersects(ray,box,thickness);
    int outside_box=box.Outside(ray.endpoint,thickness);

    if(outside_box && !intersect_box) return false; // missed the box
    else if(outside_box){ // intersected the box from the outside
        VECTOR<T,3> point=ray.Point(ray.t_max); // intersection point with the box
        point=box.Thickened(-4*thickness).Clamp(point); // moves the point inside the box
        T phi_temp;
        region_end=levelset_multiple.Inside_Region(point,phi_temp);
        region_start=-1;    // -1 means outside of levelset_multiple_object
        return true;}
    else if(!intersect_box){t_start=0;t_end=ray.t_max;} // intersected some object inside the box
    else{ // intersects the box from inside
        t_start=0;t_end=ray.t_max;
        exit_intersection=true;exit_t_max=ray.t_max;exit_aggregate=ray.aggregate_id; // save for exiting rays
        ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;} // box intersection doesn't count    

    // start marching
    T t1=t_start+thickness;
    T phi1;region_start=levelset_multiple.Inside_Region(ray.Point(t1),phi1);
    T t2=t1+Integration_Step(phi1);
    // march through the line segment
    while(t2 <= t_end){
        T phi2;region_end=levelset_multiple.Inside_Region(ray.Point(t2),phi2);
        if(region_start != region_end){
            LEVELSET_MULTIPLE_ON_A_RAY<T_LEVELSET_MULTIPLE> levelset_multiple_on_a_ray(levelset_multiple,ray,region_start);
            ITERATIVE_SOLVER<T> iterative_solver;
            iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
            ray.semi_infinite=false;
            ray.t_max=iterative_solver.Bisection_Secant_Root(levelset_multiple_on_a_ray,t1,t2);
            ray.aggregate_id=-1;
            return true;}
        else{
            t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}
    // check the last piece of the line segment
    t2=t_end;
    T phi2;region_end=levelset_multiple.Inside_Region(ray.Point(t_end),phi2);
    if(region_start != region_end){
        LEVELSET_MULTIPLE_ON_A_RAY<T_LEVELSET_MULTIPLE> levelset_multiple_on_a_ray(levelset_multiple,ray,region_start);
        ITERATIVE_SOLVER<T> iterative_solver;
        iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
        ray.semi_infinite=false;
        ray.t_max=iterative_solver.Bisection_Secant_Root(levelset_multiple_on_a_ray,t1,t2);
        ray.aggregate_id=-1;
        return true;}
    else if(exit_intersection){
        region_end=-1;ray.semi_infinite=false;ray.t_max=exit_t_max;ray.aggregate_id=exit_aggregate;
        return true;}
    else  return false;// exiting ray
}
//#####################################################################
template class RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<GRID<VECTOR<float,3> > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<GRID<VECTOR<double,3> > > >;
#endif
