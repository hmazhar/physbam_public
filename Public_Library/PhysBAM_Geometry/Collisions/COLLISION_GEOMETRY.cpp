//#####################################################################
// Copyright 2006-2008, Elliot English, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_PARTICLE_STATE.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COLLISION_GEOMETRY<TV>::
COLLISION_GEOMETRY()
    :active(true),add_to_spatial_partition(true),collision_geometry_id(0),impulse_accumulator(0),collision_geometries_for_rasterization(0)
{
    Set_Collision_Thickness();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COLLISION_GEOMETRY<TV>::
~COLLISION_GEOMETRY()
{
    delete impulse_accumulator;
}
//#####################################################################
// Function Get_Body_Penetration
//#####################################################################
template<class TV> bool COLLISION_GEOMETRY<TV>::
Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,T& hit_time,int& simplex_id,T& start_phi,T& end_phi,TV& end_body_normal,TV& body_velocity) const
{
    T temp_hit_time;TV closest_point;T_WEIGHTS weights;bool got_closest_simplex=false;
    if(Earliest_Simplex_Crossover(start_X,end_X,dt,temp_hit_time,weights,simplex_id)){
        closest_point=Simplex_World_Space_Point_From_Barycentric_Coordinates(simplex_id,weights);
        got_closest_simplex=true;hit_time=temp_hit_time;}
    else{
        // We should probably use simplex positions but this would be more expensive since we need an initialized hierarchy
        T distance;
        closest_point=Simplex_Closest_Point_On_Boundary(end_X,contour_value,(T)0,&simplex_id,&distance);
        if(distance<contour_value){got_closest_simplex=true;hit_time=dt;}}
    if(got_closest_simplex){
        T_SIMPLEX start_simplex=World_Space_Simplex(simplex_id,false);
        T_SIMPLEX end_simplex=World_Space_Simplex(simplex_id,true);
        end_body_normal=end_simplex.Normal();
        end_phi=TV::Dot_Product(end_body_normal,end_X-end_simplex.x1);
        start_phi=TV::Dot_Product(start_simplex.Normal(),start_X-start_simplex.x1);
        if(start_phi<0){end_body_normal*=-1;end_phi*=-1;start_phi*=-1;}
        body_velocity=Pointwise_Object_Velocity(simplex_id,closest_point);
        return true;}
    else return false;
}
//#####################################################################
// Function Push_Out_Point
//#####################################################################
template<class TV> bool COLLISION_GEOMETRY<TV>::
Push_Out_Point(TV& X,const T collision_distance,T& distance) const
{
    int simplex_id;
    Simplex_Closest_Point_On_Boundary(X,collision_distance,(T)0,&simplex_id,&distance);
    if(distance<collision_distance){
        T_SIMPLEX simplex=World_Space_Simplex(simplex_id,false);
        TV normal=simplex.Normal();
        T old_distance=TV::Dot_Product(X-simplex.x1,normal);
        X+=(sign(old_distance)*collision_distance-old_distance)*normal;
        return true;}
    return false;
}
//#####################################################################
// Undefined Functions
//#####################################################################
template<class TV> TV COLLISION_GEOMETRY<TV>::Pointwise_Object_Velocity(const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV COLLISION_GEOMETRY<TV>::Pointwise_Object_Velocity(const int aggregate_id,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Update_Bounding_Box() {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> const RANGE<TV>& COLLISION_GEOMETRY<TV>::Axis_Aligned_Bounding_Box() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Has_Volumetric_Geometry() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

// interface independent
template<class TV> TV COLLISION_GEOMETRY<TV>::Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2,const int max_iterations,int* simplex_id,T* distance) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

// levelset interface
template<class TV> typename TV::SCALAR COLLISION_GEOMETRY<TV>::Implicit_Geometry_Extended_Value(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV COLLISION_GEOMETRY<TV>::Implicit_Geometry_Normal(const TV& location,const int aggregate) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
// TODO: remove need for location_particle_index parameter
template<class TV> TV COLLISION_GEOMETRY<TV>::Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV COLLISION_GEOMETRY<TV>::Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV COLLISION_GEOMETRY<TV>::Implicit_Geometry_Closest_Point_On_Boundary(const TV& location,const T tolerance,const int max_iterations,T* distance) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Implicit_Geometry_Lazy_Inside_Extended_Levelset(const TV& location,T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

// simplex interface
template<class TV> TV COLLISION_GEOMETRY<TV>::Simplex_World_Space_Point_From_Barycentric_Coordinates(const int simplex_id,const T_WEIGHTS& weights) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> int COLLISION_GEOMETRY<TV>::Number_Of_Simplices() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Simplex_Intersection(RAY<TV>& ray) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Inside_Any_Simplex(const TV& location,int& simplex_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Inside(const TV& location,const T thickness_over_two) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV COLLISION_GEOMETRY<TV>::Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2,int* simplex_id,T* distance) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV COLLISION_GEOMETRY<TV>::Pointwise_Object_Pseudo_Velocity(const int aggregate_id,const TV& X,const int state1,const int state2) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename COLLISION_GEOMETRY<TV>::T_SIMPLEX COLLISION_GEOMETRY<TV>::World_Space_Simplex(const int simplex_id,const bool use_saved_state) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id,
    POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool COLLISION_GEOMETRY<TV>::Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Update_Intersection_Acceleration_Structures(const bool use_swept_triangle_hierarchy,const int state1,const int state2){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Save_State(const int state_index,const T time){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Restore_State(const int state_index){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Average_States(const int state1,const int state2,const int result_state,const T interpolation_distance){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Delete_State(const int state_index){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class TV> void COLLISION_GEOMETRY<TV>::Read_State(TYPED_ISTREAM& input,const int state_index){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void COLLISION_GEOMETRY<TV>::Write_State(TYPED_OSTREAM& output,const int state_index) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
#endif
//#####################################################################
template class COLLISION_GEOMETRY<VECTOR<float,1> >;
template class COLLISION_GEOMETRY<VECTOR<float,2> >;
template class COLLISION_GEOMETRY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COLLISION_GEOMETRY<VECTOR<double,1> >;
template class COLLISION_GEOMETRY<VECTOR<double,2> >;
template class COLLISION_GEOMETRY<VECTOR<double,3> >;
#endif
