//#####################################################################
// Copyright 2007, Andrew Selle, Justin Solomon.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_SEGMENTED_CURVE
//##################################################################### 
#ifndef __RENDERING_SEGMENTED_CURVE__
#define __RENDERING_SEGMENTED_CURVE__

#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_CYLINDER_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_SEGMENTED_CURVE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;
    using RENDERING_OBJECT<T>::World_Space_Bounding_Box;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    SEGMENTED_CURVE<TV>& segmented_curve;
    bool add_segments_to_acceleration_structure;
    T thickness;

    RENDERING_SEGMENTED_CURVE(SEGMENTED_CURVE<TV>& segmented_curve_input, T thickness_input=1e-4)
        :segmented_curve(segmented_curve_input),add_segments_to_acceleration_structure(true),thickness(thickness_input)
    {
        segmented_curve.Update_Bounding_Box();
        segmented_curve.bounding_box->Change_Size(thickness); // thickness accounts for cylinder size
        if(!segmented_curve.hierarchy) segmented_curve.Initialize_Hierarchy(true); // TODO: when spatial partition, do we need hierachy?
    }

    virtual ~RENDERING_SEGMENTED_CURVE()
    {}

    bool Intersection(RAY<TV> &ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    bool intersected=false;
    for(int t=1;t<=segmented_curve.mesh.elements.m;t++){ // TODO: make this more efficient using hierarchy.
        const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(t);
        if(INTERSECTION::Intersects(ray,CYLINDER<T>(segmented_curve.particles.X(nodes[1]),segmented_curve.particles.X(nodes[2]),thickness))) {
            ray.aggregate_id=t;intersected=true;}}
    return intersected;}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    return World_Space_Vector(CYLINDER<T>(segmented_curve.particles.X(nodes[1]),segmented_curve.particles.X(nodes[2]),thickness).Normal(Object_Space_Point(location),1));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return false;} // segmented curves have no inside

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return true;}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {T distance;segmented_curve.Closest_Point_On_Curve(Object_Space_Point(location),0,NULL,&distance);
    return (abs(distance-thickness) < small_number);}

    bool Intersection(RAY<TV>& ray,const int aggregate) const PHYSBAM_OVERRIDE
    {if(!add_segments_to_acceleration_structure) return Intersection(ray);
    const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    if(INTERSECTION::Intersects(ray,CYLINDER<T>(segmented_curve.particles.X(nodes[1]),segmented_curve.particles.X(nodes[2]),thickness))){
        ray.aggregate_id=aggregate;return true;}
    return false;}

    void Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const PHYSBAM_OVERRIDE
    {if(add_segments_to_acceleration_structure){
        for(int t=1;t<=segmented_curve.mesh.elements.m;t++){
            const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(t);
            RANGE<TV> bounding_box=RANGE<TV>::Bounding_Box(World_Space_Point(segmented_curve.particles.X(nodes[1])),World_Space_Point(segmented_curve.particles.X(nodes[2])));
            bounding_box.Change_Size(thickness);
            primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(bounding_box,this,t));}}
    else{
        primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(Object_Space_Bounding_Box(),this,1));}}

    bool Has_Bounding_Box() const PHYSBAM_OVERRIDE
    {return true;}

    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {if(!segmented_curve.bounding_box){
        segmented_curve.Update_Bounding_Box();
        segmented_curve.bounding_box->Change_Size(thickness);} // thickness accounts for cylinder size
    return *segmented_curve.bounding_box;}

    void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const PHYSBAM_OVERRIDE
    {const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    object_tangent=segmented_curve.particles.X(nodes[1])-segmented_curve.particles.X(nodes[2]).Normalized();
    object_bitangent=TV::Cross_Product(object_tangent,object_space_normal).Normalized();}

    void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const PHYSBAM_OVERRIDE
    {const VECTOR<int,2>& nodes=segmented_curve.mesh.elements(aggregate);
    world_tangent=World_Space_Vector(segmented_curve.particles.X(nodes[1])-segmented_curve.particles.X(nodes[2])).Normalized();
    world_bitangent=TV::Cross_Product(world_tangent,world_space_normal).Normalized();}

//#####################################################################
};
}
 
#endif
