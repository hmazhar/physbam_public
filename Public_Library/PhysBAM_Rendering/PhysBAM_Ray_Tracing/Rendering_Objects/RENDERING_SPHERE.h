//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_SPHERE
//##################################################################### 
#ifndef __RENDERING_SPHERE__
#define __RENDERING_SPHERE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SPHERE_INTERSECTION.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_SPHERE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::World_Space_Bounding_Box;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    SPHERE<TV> sphere;

    RENDERING_SPHERE()
    {}

    RENDERING_SPHERE(const TV& center_input,const T radius_input)
        :sphere(center_input,radius_input)
    {}

    virtual ~RENDERING_SPHERE()
    {}

    void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const PHYSBAM_OVERRIDE
    {TV p=(object_space_point-sphere.center)/sphere.radius;s=((T)pi+atan2(p.z,p.x))/T(2*pi);t=acos(p.y)/(T)pi;}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,sphere,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    bool Intersection(RAY<TV>& ray,const int aggregate) const PHYSBAM_OVERRIDE
    {return Intersection(ray);}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(sphere.Normal(Object_Space_Point(location)));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return sphere.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return sphere.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return sphere.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const PHYSBAM_OVERRIDE
    {return World_Space_Point(sphere.Surface(Object_Space_Point(location)));}

    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {return sphere.Signed_Distance(Object_Space_Point(location));}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {return TESSELLATION::Generate_Triangles(sphere);}

    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {TV radius_vector(sphere.radius,sphere.radius,sphere.radius);return RANGE<TV>(sphere.center-radius_vector,sphere.center+radius_vector);}

    void Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const PHYSBAM_OVERRIDE
    {primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(World_Space_Bounding_Box(),this,0));}

    void Get_Object_Space_Tangent_And_Bitangent(const VECTOR<T,3>& object_space_point,const VECTOR<T,3>& object_space_normal,const int aggregate,VECTOR<T,3>& object_tangent,VECTOR<T,3>& object_bitangent) const PHYSBAM_OVERRIDE
    {T x=object_space_point.x-sphere.center.x;T y=object_space_point.y-sphere.center.y;T z=object_space_point.z-sphere.center.z;
    T theta=asin(z/sphere.radius);T phi=asin(y/sqrt(x*x+y*y));
    if ((T)cos(phi)*(T)sqrt(sqr(x)+sqr(y))-x > (T)1e-3) phi=(T)pi-phi;
    object_tangent.Set(-sphere.radius*sin(phi)*cos(theta),sphere.radius*cos(theta)*cos(phi),0);object_tangent.Normalize();
    object_bitangent=VECTOR<T,3>::Cross_Product(object_tangent,object_space_normal);object_bitangent.Normalize();}

    void Get_World_Space_Tangent_And_Bitangent(const VECTOR<T,3>& world_space_point,const VECTOR<T,3>& world_space_normal,const int aggregate,VECTOR<T,3>& world_tangent,VECTOR<T,3>& world_bitangent) const PHYSBAM_OVERRIDE
    {TV object_space_point=Object_Space_Point(world_space_point);
    TV object_space_normal=Object_Space_Vector(world_space_normal);
    Get_Object_Space_Tangent_And_Bitangent(object_space_point,object_space_normal,aggregate,world_tangent,world_bitangent);
    world_tangent=World_Space_Vector(world_tangent);world_bitangent=World_Space_Vector(world_bitangent);}
    
//#####################################################################
};    
}
#endif
