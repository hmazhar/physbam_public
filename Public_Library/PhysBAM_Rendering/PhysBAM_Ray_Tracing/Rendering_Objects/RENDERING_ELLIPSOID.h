//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_ELLIPSOID
//#####################################################################
#ifndef __RENDERING_ELLIPSOID__
#define __RENDERING_ELLIPSOID__

#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_ELLIPSOID:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;

    ELLIPSOID<T> ellipsoid;

    RENDERING_ELLIPSOID()
    {}

    RENDERING_ELLIPSOID(const TV& center_input,const TV& radius_input,const ROTATION<TV >& orientation_input)
        :ellipsoid(center_input,radius_input,orientation_input)
    {}

    virtual ~RENDERING_ELLIPSOID()
    {}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(ellipsoid.Intersection(object_space_ray,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(ellipsoid.Normal(Object_Space_Point(location)));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return ellipsoid.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return ellipsoid.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return ellipsoid.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const PHYSBAM_OVERRIDE
    {return World_Space_Point(ellipsoid.Approximate_Surface(Object_Space_Point(location),small_number));}

    // PHYSBAM_WARNING - THIS CALLS APPROXIMATE SIGNED DISTANCE!!!!
    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {return ellipsoid.Approximate_Signed_Distance(Object_Space_Point(location));}

//#####################################################################
};
}
#endif

