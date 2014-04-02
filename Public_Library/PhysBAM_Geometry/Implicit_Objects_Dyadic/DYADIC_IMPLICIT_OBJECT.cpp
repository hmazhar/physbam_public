//#####################################################################
// Copyright 2003-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BLOCK_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/DYADIC_IMPLICIT_OBJECT_ON_A_RAY.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
using namespace PhysBAM;
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Dyadic_Implicit_Object(){
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<DYADIC_IMPLICIT_OBJECT<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<DYADIC_IMPLICIT_OBJECT<VECTOR<float,3> > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<DYADIC_IMPLICIT_OBJECT<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<DYADIC_IMPLICIT_OBJECT<VECTOR<double,3> > >();
#endif
    return true;
}
static bool registered=Register_Dyadic_Implicit_Object();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DYADIC_IMPLICIT_OBJECT<TV>::
DYADIC_IMPLICIT_OBJECT(T_GRID& grid_input,ARRAY<T>& phi_input,ARRAY<T>* phi_nodes_input)
    :levelset(grid_input,phi_input),phi_nodes(phi_nodes_input),need_destroy_data(false)
{
    PHYSBAM_ASSERT(registered);
    Update_Box();Update_Minimum_Cell_Size(levelset.grid.maximum_depth);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DYADIC_IMPLICIT_OBJECT<TV>::
~DYADIC_IMPLICIT_OBJECT()
{
    if(need_destroy_data){delete &levelset.grid;delete &levelset.phi;}
    delete phi_nodes;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> DYADIC_IMPLICIT_OBJECT<TV>* DYADIC_IMPLICIT_OBJECT<TV>::
Create()
{
    DYADIC_IMPLICIT_OBJECT* implicit_object=new DYADIC_IMPLICIT_OBJECT(*(new T_GRID),*(new ARRAY<T>));
    implicit_object->need_destroy_data=true;return implicit_object;
}
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
template<class TV> bool DYADIC_IMPLICIT_OBJECT<TV>::
Intersection(RAY<TV>& ray,const T thickness,const ARRAY<T>& phi_nodes,const bool verbose) const
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;

    T t_start,t_end; 
    bool exit_intersection=false;int exit_aggregate=0;T exit_t_max=0;
    int intersect_box=INTERSECTION::Intersects(ray,box,thickness);
    int outside_box=box.Outside(ray.endpoint,thickness);

    if(outside_box && !intersect_box) return false; // missed the box
    else if(outside_box){ // intersected the box from the outside
        TV point=ray.Point(ray.t_max); // intersection point with the box
        if((*this)(point) <= 0) return true; // test the levelset right at the box
        point=box.Thickened(-4*thickness).Clamp(point); // moves the point inside the box
        if((*this)(point) <= 0) return true; // level set is on the edge of the box
        else{ // level set is not on the edge of the box
            t_start=ray.t_max; // intersection with the box
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate; // box intersection doesn't count
            RAY<TV> new_ray(point,ray.direction);INTERSECTION::Intersects(new_ray,box,thickness);
            if(ray.semi_infinite) t_end=t_start+new_ray.t_max;
            else t_end=min(t_start+new_ray.t_max,ray.t_max);}}
    else if(!intersect_box){t_start=0;t_end=ray.t_max;} // intersected some object inside the box
    else{ // intersects the box from inside
        t_start=0;t_end=ray.t_max;
        exit_intersection=true;exit_t_max=ray.t_max;exit_aggregate=ray.aggregate_id; // save for exiting rays
        ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;} // box intersection doesn't count

    // set up marching
    DYADIC_IMPLICIT_OBJECT_ON_A_RAY<T_GRID> implicit_surface_on_a_ray(*this,ray);
    ITERATIVE_SOLVER<T> iterative_solver;iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
    // start marching
    T t1=t_start+thickness;
    
    TV p1(ray.Point(t1));T_CELL* current_cell=levelset.grid.Leaf_Cell(p1);
    T phi1=levelset.Phi_From_Close_Cell(current_cell,p1,&phi_nodes),t2=t1+Integration_Step(phi1);
    // march through the line segment
    while(t2 <= t_end){
        TV p2(ray.Point(t2));
        if(!levelset.grid.Inside_Cell(current_cell,p2)){
            if(phi1<current_cell->DX().x*20) // 20 is a hack to avoid looking for neighbors
                current_cell=levelset.grid.Leaf_Cell_By_Neighbor_Path(current_cell,p2);
            else current_cell=levelset.grid.Leaf_Cell(p2);}
        else if(current_cell->Has_Children()) current_cell=levelset.grid.Leaf_Cell(p2);
        T phi2=levelset.Phi_From_Close_Cell(current_cell,p2,&phi_nodes);
        if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
            T_CELL* base_cell=levelset.grid.Base_Cell_By_Neighbor_Path(current_cell,p1);assert(base_cell);
            implicit_surface_on_a_ray.phi_nodes=&phi_nodes;
            implicit_surface_on_a_ray.last_block=new BLOCK_DYADIC<T_GRID>(levelset.grid,base_cell); // optimization
            ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;
            if(implicit_surface_on_a_ray.last_block){delete implicit_surface_on_a_ray.last_block;implicit_surface_on_a_ray.last_block=0;}
            return true;}
        else{t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}} 
    // check the last piece of the line segment
    t2=t_end;T phi2=(*this)(ray.Point(t_end));
    if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;return true;}    
    else if(exit_intersection && phi2 <= 0){ray.semi_infinite=false;ray.t_max=exit_t_max; ray.aggregate_id=exit_aggregate; return true;} // exiting ray
    else return false;
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template DYADIC_IMPLICIT_OBJECT<VECTOR<T,d> >::DYADIC_IMPLICIT_OBJECT(DYADIC_GRID_POLICY<VECTOR<T,d> >::DYADIC_GRID& grid_input,ARRAY<T>& phi_input,ARRAY<T>* phi_nodes_input); \
    template DYADIC_IMPLICIT_OBJECT<VECTOR<T,d> >* DYADIC_IMPLICIT_OBJECT<VECTOR<T,d> >::Create();
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
template bool DYADIC_IMPLICIT_OBJECT<VECTOR<float,3> >::Intersection(RAY<VECTOR<float,3> >&,const float,const ARRAY<float>&,const bool) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
template bool DYADIC_IMPLICIT_OBJECT<VECTOR<double,3> >::Intersection(RAY<VECTOR<double,3> >&,const double,const ARRAY<double>&,const bool) const;
#endif
#endif
