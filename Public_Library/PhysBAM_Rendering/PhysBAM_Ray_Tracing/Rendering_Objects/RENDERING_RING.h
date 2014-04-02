//#####################################################################
// Copyright 2003, Eran Guendelman, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_RING
//#####################################################################
#ifndef __RENDERING_RING__
#define __RENDERING_RING__

#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_RING:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;

    RING<T> ring;

    RENDERING_RING()
    {}

    RENDERING_RING(const TV& point1_input,const TV& point2_input,const T outer_radius_input,const T inner_radius_input)
        :ring(point1_input,point2_input,outer_radius_input,inner_radius_input)
    {}

    virtual ~RENDERING_RING()
    {}

    bool Intersection(RAY<TV>& ray) const    PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(ring.Intersection(object_space_ray,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(ring.Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return ring.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return ring.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return ring.Boundary(Object_Space_Point(location),small_number);}

    bool Has_Bounding_Box() const  PHYSBAM_OVERRIDE
    {return true;}
    
    // This is wrong in general!!!  Only works if the two ring points are aligned along the y axis and have x=z=0
    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {return RANGE<TV>(-ring.outer_radius,ring.outer_radius,min(ring.plane1.x1.y,ring.plane2.x1.y),max(ring.plane1.x1.y,ring.plane2.x1.y),-ring.outer_radius,ring.outer_radius);}

//#####################################################################
};   
}
#endif

