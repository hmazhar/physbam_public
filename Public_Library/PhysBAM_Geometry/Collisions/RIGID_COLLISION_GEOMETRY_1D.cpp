//#####################################################################
// Copyright 2007, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
RIGID_COLLISION_GEOMETRY(RIGID_GEOMETRY<TV>& rigid_geometry_input)
    :RIGID_COLLISION_GEOMETRY_BASE<TV>(rigid_geometry_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
~RIGID_COLLISION_GEOMETRY()
{}
//##################################################################### 
// Function Earliest_Simplex_Crossover
//##################################################################### 
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
Earliest_Simplex_Crossover(const VECTOR<T,1>& start_X,const VECTOR<T,1>& end_X,const T dt,T& hit_time,ONE& weight,int& segment_id) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    T min_time=FLT_MAX;bool collision=false;T current_hit_time,relative_speed;ONE current_weight;
    TV normal;
    for(int segment_number=1;segment_number<=rigid_geometry.simplicial_object->mesh.elements.m;segment_number++){
        POINT_SIMPLEX_1D<T> initial_segment=World_Space_Simplex(segment_number),final_segment=World_Space_Simplex(segment_number,saved_states(1).x);
        POINT_SIMPLEX_COLLISION_TYPE collision_type=
            CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Point_Collision(initial_segment,final_segment,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,normal,current_weight,relative_speed);
        if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time<min_time){
            min_time=hit_time=current_hit_time;weight=current_weight;segment_id=segment_number;collision=true;}}
    return collision;
}
//##################################################################### 
// Function Latest_Simplex_Crossover
//##################################################################### 
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
Latest_Simplex_Crossover(const VECTOR<T,1>& start_X,const VECTOR<T,1>& end_X,const T dt,T& hit_time,ONE& weight,int& segment_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    returned_collision_type=POINT_SIMPLEX_NO_COLLISION;
    T max_time=-FLT_MAX;bool collision=false;T current_hit_time,relative_speed;ONE current_weight;
    TV normal;
    for (int segment_number=1;segment_number<=rigid_geometry.simplicial_object->mesh.elements.m;segment_number++){
        POINT_SIMPLEX_1D<T> initial_segment=World_Space_Simplex(segment_number),final_segment=World_Space_Simplex(segment_number,saved_states(1).x);
        POINT_SIMPLEX_COLLISION_TYPE collision_type=
            CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Point_Collision(initial_segment,final_segment,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,normal,current_weight,relative_speed);
        if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time>max_time){
            max_time=hit_time=current_hit_time;weight=current_weight;segment_id=segment_number;collision=true;returned_collision_type=collision_type;}}
    return collision;
}
//#####################################################################
// Function Any_Simplex_Crossover
//#####################################################################
template<class T> bool RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
Any_Simplex_Crossover(const VECTOR<T,1>& start_X,const VECTOR<T,1>& end_X,const T dt) const
{
    const T collision_thickness_over_two=(T).5*collision_thickness;
    T current_hit_time,relative_speed;ONE current_weight;
    TV normal;
    for(int segment_number=1;segment_number<=rigid_geometry.simplicial_object->mesh.elements.m;segment_number++){
        POINT_SIMPLEX_1D<T> initial_segment=World_Space_Simplex(segment_number),final_segment=World_Space_Simplex(segment_number,saved_states(1).x);
        POINT_SIMPLEX_COLLISION_TYPE collision_type=
            CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Point_Collision(initial_segment,final_segment,start_X,end_X,dt,collision_thickness_over_two,current_hit_time,normal,current_weight,relative_speed);
        if(collision_type!=POINT_SIMPLEX_NO_COLLISION) return true;}
    return false;
}
//##################################################################### 
// Function Get_Simplex_Bounding_Boxes
//##################################################################### 
template<class T> void RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const
{
    if(!rigid_geometry.simplicial_object->point_simplex_list) rigid_geometry.simplicial_object->Update_Point_Simplex_List();
    for(int t=1;t<=rigid_geometry.simplicial_object->mesh.elements.m;t++){
        RANGE<TV> box=rigid_geometry.World_Space_Simplex_Bounding_Box(t);
        if(with_body_motion) box.Enlarge_To_Include_Box(rigid_geometry.World_Space_Simplex_Bounding_Box(t,saved_states(1).x));
        box.Change_Size(extra_thickness+body_thickness_factor*collision_thickness);
        bounding_boxes.Append(box);}
}
//##################################################################### 
// Function World_Space_Simplex
//##################################################################### 
template<class T> POINT_SIMPLEX_1D<T> RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
World_Space_Simplex(const int segment_id,const bool use_saved_state) const
{
    if(use_saved_state){return World_Space_Simplex(segment_id,saved_states(1).x);}
    return rigid_geometry.World_Space_Simplex(segment_id);
}
//##################################################################### 
// Function World_Space_Simplex
//##################################################################### 
template<class T> POINT_SIMPLEX_1D<T> RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
World_Space_Simplex(const int segment_id,const FRAME<TV>& state) const
{
    POINT_SIMPLEX_1D<T>& object_space_point_simplex=(*rigid_geometry.simplicial_object->point_simplex_list)(segment_id);
    return POINT_SIMPLEX_1D<T>(state*object_space_point_simplex.x1,object_space_point_simplex.direction);
}
//#####################################################################
// Function Update_Intersection_Acceleration_Structures
//#####################################################################
template<class T> void RIGID_COLLISION_GEOMETRY<VECTOR<T,1> >::
Update_Intersection_Acceleration_Structures(const bool use_swept_triangle_hierarchy,const int state1,const int state2)
{}
//##################################################################### 
template class RIGID_COLLISION_GEOMETRY<VECTOR<float,1> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_COLLISION_GEOMETRY<VECTOR<double,1> >;
#endif
