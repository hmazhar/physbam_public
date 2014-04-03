//#####################################################################
// Copyright 2002, 2003, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_TETRAHEDRON
//#####################################################################
#ifndef __RENDERING_TETRAHEDRON__
#define __RENDERING_TETRAHEDRON__

#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_TETRAHEDRON:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;

    TETRAHEDRON<T> tetrahedron;

    RENDERING_TETRAHEDRON()
    {}

    RENDERING_TETRAHEDRON(TV& x1_input,TV& x2_input,TV& x3_input,TV& x4_input)
        :tetrahedron(x1_input,x2_input,x3_input,x4_input)
    {}

    virtual ~RENDERING_TETRAHEDRON()
    {}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(tetrahedron.Intersection(object_space_ray,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(tetrahedron.Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return tetrahedron.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return tetrahedron.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return tetrahedron.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const PHYSBAM_OVERRIDE
    {return World_Space_Point(tetrahedron.Surface(Object_Space_Point(location)));}

//#####################################################################
};   
}
#endif

