//#####################################################################
// Copyright 2006-2008, Elliot English, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#endif
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_POINT_SIMPLICES_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
//using namespace PhysBAM;
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_COLLISION_GEOMETRY_BASE<TV>::
RIGID_COLLISION_GEOMETRY_BASE(RIGID_GEOMETRY<TV>& rigid_geometry_input)
    :rigid_geometry(rigid_geometry_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_COLLISION_GEOMETRY_BASE<TV>::
~RIGID_COLLISION_GEOMETRY_BASE()
{
}
//#####################################################################
// Function Pointwise_Object_Velocity
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Pointwise_Object_Velocity(const TV& X) const
{
    return rigid_geometry.Pointwise_Object_Velocity(X);
}
//#####################################################################
// Function Pointwise_Object_Velocity
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Pointwise_Object_Velocity(const int simplex_id,const TV& X) const // extra simplex_id is not used, but for a virtual function in COLLISION GEOMETRY
{
    return Pointwise_Object_Velocity(X);
}
//#####################################################################
// Function Pointwise_Object_Velocity_At_Particle
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Pointwise_Object_Velocity_At_Particle(const TV& X,const int particle_index) const
{
    return rigid_geometry.Pointwise_Object_Velocity_At_Particle(X,particle_index);
}
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
template<class TV> void RIGID_COLLISION_GEOMETRY_BASE<TV>::
Update_Bounding_Box()
{
    rigid_geometry.Update_Bounding_Box();
}
//#####################################################################
// Function Axis_Aligned_Bounding_Box
//#####################################################################
template<class TV> const RANGE<TV>& RIGID_COLLISION_GEOMETRY_BASE<TV>::
Axis_Aligned_Bounding_Box() const
{
    return rigid_geometry.Axis_Aligned_Bounding_Box();
}
//#####################################################################
// Function Has_Volumetric_Geometry
//#####################################################################
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Has_Volumetric_Geometry() const
{
    return rigid_geometry.template Find_Structure<typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT*>()!=0;
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2,const int max_iterations,int* simplex_id,T* distance) const
{
    if(rigid_geometry.implicit_object && rigid_geometry.simplicial_object){
        TV normal=Implicit_Geometry_Normal(location);
        RAY<TV> ray(location,normal);

        if(INTERSECTION::Intersects(ray,*rigid_geometry.simplicial_object,thickness_over_2)){
            if(simplex_id) *simplex_id=ray.aggregate_id;
            if(distance) *distance=ray.t_max;
            return ray.Point(ray.t_max);}
        else 
            return Simplex_Closest_Point_On_Boundary(location,max_distance,collision_thickness,simplex_id,distance);
    }
    else if(rigid_geometry.implicit_object)
        return Implicit_Geometry_Closest_Point_On_Boundary(location,thickness_over_2,max_iterations,distance);
    return Simplex_Closest_Point_On_Boundary(location,max_distance,collision_thickness,simplex_id,distance);
}
//#####################################################################
// Function Implicit_Geometry_Extended_Value
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Extended_Value(const TV& location) const
{
    return rigid_geometry.Implicit_Geometry_Extended_Value(location);
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Normal(const TV& location,const int aggregate) const
{
    return rigid_geometry.Implicit_Geometry_Normal(location,aggregate);
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    return rigid_geometry.Implicit_Geometry_Normal(location,phi_value,aggregate,location_particle_index);
}
//#####################################################################
// Function Implicit_Geometry_Extended_Normal
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    return rigid_geometry.Implicit_Geometry_Extended_Normal(location,phi_value,aggregate,location_particle_index);
}
//#####################################################################
// Function Implicit_Geometry_Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Closest_Point_On_Boundary(const TV& location,const T tolerance,const int max_iterations,T* distance) const
{
    assert(rigid_geometry.implicit_object);return rigid_geometry.implicit_object->Closest_Point_On_Boundary(location,tolerance,max_iterations,distance);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside
//#####################################################################
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value) const
{
    return rigid_geometry.Implicit_Geometry_Lazy_Inside(location,contour_value);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value) const
{
    return rigid_geometry.Implicit_Geometry_Lazy_Inside_And_Value(location,phi,contour_value);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Implicit_Geometry_Lazy_Inside_Extended_Levelset(const TV& location,T contour_value) const
{
    return rigid_geometry.implicit_object->Lazy_Inside_Extended_Levelset(location,contour_value);
}
//#####################################################################
// Function Simplex_World_Space_Point_From_Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,1> Simplex_World_Space_Point_From_Barycentric_Coordinates_Helper(RIGID_GEOMETRY<VECTOR<T,1> >& rigid_geometry,const int point_simplex_id,const ONE&)
{
    const POINT_SIMPLEX_1D<T>& object_space_point_simplex=(*rigid_geometry.simplicial_object->point_simplex_list)(point_simplex_id);
    return rigid_geometry.World_Space_Point(object_space_point_simplex.x1);
}
template<class T> VECTOR<T,2> Simplex_World_Space_Point_From_Barycentric_Coordinates_Helper(RIGID_GEOMETRY<VECTOR<T,2> >& rigid_geometry,const int segment_id,const T& alpha)
{
    const SEGMENT_2D<T>& object_space_segment=(*rigid_geometry.simplicial_object->segment_list)(segment_id);
    return SEGMENT_2D<T>::Point_From_Barycentric_Coordinates(alpha,rigid_geometry.World_Space_Point(object_space_segment.x1),rigid_geometry.World_Space_Point(object_space_segment.x2));
}
template<class T> VECTOR<T,3> Simplex_World_Space_Point_From_Barycentric_Coordinates_Helper(RIGID_GEOMETRY<VECTOR<T,3> >& rigid_geometry,const int triangle_id,const VECTOR<T,3>& weights)
{
    return TRIANGLE_3D<T>::Point_From_Barycentric_Coordinates(weights,
        rigid_geometry.World_Space_Point(rigid_geometry.simplicial_object->particles.X(rigid_geometry.simplicial_object->mesh.elements(triangle_id)(1))),
        rigid_geometry.World_Space_Point(rigid_geometry.simplicial_object->particles.X(rigid_geometry.simplicial_object->mesh.elements(triangle_id)(2))),
        rigid_geometry.World_Space_Point(rigid_geometry.simplicial_object->particles.X(rigid_geometry.simplicial_object->mesh.elements(triangle_id)(3))));
}
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Simplex_World_Space_Point_From_Barycentric_Coordinates(const int simplex_id,const T_WEIGHTS& weights) const
{
    return Simplex_World_Space_Point_From_Barycentric_Coordinates_Helper(rigid_geometry,simplex_id,weights);
}
//#####################################################################
// Function Simplex_Intersection
//#####################################################################
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Simplex_Intersection(RAY<TV>& ray) const
{
    return rigid_geometry.Simplex_Intersection(ray,(T).5*collision_thickness);
}
//##################################################################### 
// Function Simplex_Closest_Non_Intersecting_Point
//##################################################################### 
template<class T> bool Simplex_Closest_Non_Intersecting_Point_Helper(RIGID_GEOMETRY<VECTOR<T,1> >& rigid_geometry,RAY<VECTOR<T,1> >& ray,const T collision_thickness)
{
    RAY<VECTOR<T,1> > object_space_ray=rigid_geometry.Object_Space_Ray(ray);
    if(INTERSECTION::Closest_Non_Intersecting_Point(object_space_ray,*rigid_geometry.simplicial_object,collision_thickness)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    return false;
}
template<class T> bool Simplex_Closest_Non_Intersecting_Point_Helper(RIGID_GEOMETRY<VECTOR<T,2> >& rigid_geometry,RAY<VECTOR<T,2> >& ray,const T collision_thickness)
{
    RAY<VECTOR<T,2> > object_space_ray=rigid_geometry.Object_Space_Ray(ray);
    if(INTERSECTION::Closest_Non_Intersecting_Point(object_space_ray,*rigid_geometry.simplicial_object,collision_thickness)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    return false;
}
template<class T> bool Simplex_Closest_Non_Intersecting_Point_Helper(RIGID_GEOMETRY<VECTOR<T,3> >& rigid_geometry,RAY<VECTOR<T,3> >& ray,const T collision_thickness)
{
    RAY<VECTOR<T,3> > object_space_ray=rigid_geometry.Object_Space_Ray(ray);
    if(INTERSECTION::Closest_Non_Intersecting_Point(object_space_ray,*rigid_geometry.simplicial_object,collision_thickness)){
        ray.Restore_Intersection_Information(object_space_ray);return true;}
    return false;
}
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const
{
    return Simplex_Closest_Non_Intersecting_Point_Helper(rigid_geometry,ray,collision_thickness);
}
//##################################################################### 
// Function Inside_Any_Simplex
//##################################################################### 
template<class T> bool Inside_Any_Simplex_Helper(RIGID_GEOMETRY<VECTOR<T,1> >& rigid_geometry,const VECTOR<T,1>& location,int& segment_id,const T collision_thickness)
{
    return rigid_geometry.simplicial_object->Inside_Any_Simplex(rigid_geometry.Object_Space_Point(location),segment_id,collision_thickness);
}
template<class T> bool Inside_Any_Simplex_Helper(RIGID_GEOMETRY<VECTOR<T,2> >& rigid_geometry,const VECTOR<T,2>& location,int& segment_id,const T collision_thickness)
{
    return rigid_geometry.simplicial_object->Inside_Any_Segment(rigid_geometry.Object_Space_Point(location),segment_id,collision_thickness);
}
template<class T> bool Inside_Any_Simplex_Helper(RIGID_GEOMETRY<VECTOR<T,3> >& rigid_geometry,const VECTOR<T,3>& location,int& triangle_id,const T collision_thickness)
{
    return rigid_geometry.simplicial_object->Inside_Any_Triangle(rigid_geometry.Object_Space_Point(location),triangle_id,collision_thickness);
}
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Inside_Any_Simplex(const TV& location,int& simplex_id) const
{
    return Inside_Any_Simplex_Helper(rigid_geometry,location,simplex_id,collision_thickness);
}
//##################################################################### 
// Function Inside
//##################################################################### 
template<class TV> bool RIGID_COLLISION_GEOMETRY_BASE<TV>::
Inside(const TV& location,const T thickness_over_two) const
{
    return rigid_geometry.Simplex_Inside(location,thickness_over_two);
}
//#####################################################################
// Function Simplex_Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2,int* simplex_id,T* returned_distance) const
{
    assert(rigid_geometry.simplicial_object);
    return rigid_geometry.World_Space_Point(rigid_geometry.simplicial_object->Closest_Point_On_Boundary(rigid_geometry.Object_Space_Point(location),
            max_distance,thickness_over_2,simplex_id,returned_distance));
}
//#####################################################################
// Function Pointwise_Object_Pseudo_Velocity
//#####################################################################
template<class TV> TV RIGID_COLLISION_GEOMETRY_BASE<TV>::
Pointwise_Object_Pseudo_Velocity(const int simplex_id,const TV& X,const int state1,const int state2) const
{
    TV object_space_point=rigid_geometry.Object_Space_Point(X);
    return (saved_states(state2).x*object_space_point-saved_states(state1).x*object_space_point)/(saved_states(state2).y-saved_states(state1).y);
}
//#####################################################################
// Function Number_Of_Simplices
//#####################################################################
template<class TV> int RIGID_COLLISION_GEOMETRY_BASE<TV>::
Number_Of_Simplices() const
{
    return rigid_geometry.simplicial_object->mesh.elements.m;
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class TV> void RIGID_COLLISION_GEOMETRY_BASE<TV>::
Save_State(const int state_index,const T time)
{
    if(saved_states.m<state_index) saved_states.Resize(state_index);
    saved_states(state_index).x=rigid_geometry.Frame();
    saved_states(state_index).y=time;
}
//#####################################################################
// Function Restore_State
//#####################################################################
template<class TV> void RIGID_COLLISION_GEOMETRY_BASE<TV>::
Restore_State(const int state_index)
{
    rigid_geometry.Set_Frame(saved_states(state_index).x);
}
//#####################################################################
// Function Average_States
//#####################################################################
template<class TV> void RIGID_COLLISION_GEOMETRY_BASE<TV>::
Average_States(const int state1,const int state2,const int result_state,const T interpolation_distance)
{
    if(saved_states.m<result_state) saved_states.Resize(result_state);
    saved_states(result_state).x.t=((T)1-interpolation_distance)*saved_states(state1).x.t+interpolation_distance*saved_states(state2).x.t;
    saved_states(result_state).x.r=ROTATION<TV>::Spherical_Linear_Interpolation(saved_states(state1).x.r,saved_states(state2).x.r,interpolation_distance);
    saved_states(result_state).y=((T)1-interpolation_distance)*saved_states(state1).y+interpolation_distance*saved_states(state2).y;
}
//#####################################################################
// Function Delete_State
//#####################################################################
template<class TV> void RIGID_COLLISION_GEOMETRY_BASE<TV>::
Delete_State(const int state_index)
{
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
//#####################################################################
// Function Read_State
//#####################################################################
template<class TV> void RIGID_COLLISION_GEOMETRY_BASE<TV>::
Read_State(TYPED_ISTREAM& input,const int state_index)
{
    if(saved_states.m<state_index) saved_states.Resize(state_index);
    Read_Binary(input,saved_states(state_index));
}
//#####################################################################
// Function Write_State
//#####################################################################
template<class TV> void RIGID_COLLISION_GEOMETRY_BASE<TV>::
Write_State(TYPED_OSTREAM& output,const int state_index) const
{
    Write_Binary(output,saved_states(state_index));
}
#endif
//#####################################################################
template class RIGID_COLLISION_GEOMETRY_BASE<VECTOR<float,1> >;
template class RIGID_COLLISION_GEOMETRY_BASE<VECTOR<float,2> >;
template class RIGID_COLLISION_GEOMETRY_BASE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_COLLISION_GEOMETRY_BASE<VECTOR<double,1> >;
template class RIGID_COLLISION_GEOMETRY_BASE<VECTOR<double,2> >;
template class RIGID_COLLISION_GEOMETRY_BASE<VECTOR<double,3> >;
#endif
}
