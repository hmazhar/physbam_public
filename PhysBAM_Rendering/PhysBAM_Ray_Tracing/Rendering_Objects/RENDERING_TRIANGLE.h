//#####################################################################
// Copyright 2002, 2003, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_TRIANGLE
//#####################################################################
#ifndef __RENDERING_TRIANGLE__
#define __RENDERING_TRIANGLE__

#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_PLANE.h>
namespace PhysBAM{

template<class T>
class RENDERING_TRIANGLE:public RENDERING_PLANE<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;
    using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    TRIANGLE_3D<T> triangle;

    RENDERING_TRIANGLE()
    {}

    RENDERING_TRIANGLE(const TV& x1_input,const TV& x2_input,const TV& x3_input)
        :RENDERING_PLANE<T>(x1_input,x2_input,x3_input),triangle(x1_input,x2_input,x3_input)
    {}

    virtual ~RENDERING_TRIANGLE()
    {}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,triangle,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

//#####################################################################
};   
}
#endif

