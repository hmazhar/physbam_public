//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_BLEND_TRIANGULATED_SURFACE_SHADER__
#define __RENDERING_BLEND_TRIANGULATED_SURFACE_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BLEND_TRIANGULATED_SURFACE_SHADER:public RENDERING_BLEND_SHADER<T>
{
public:
    const ARRAY<T>& field;
    T low_value,high_value;
    
    RENDERING_BLEND_TRIANGULATED_SURFACE_SHADER(const ARRAY<T>& field_input,const T low_value_input,const T high_value_input,const MATERIAL_SHADER<T>& shader1_input,
        const MATERIAL_SHADER<T>& shader2_input,const bool direct_shading_only,RENDER_WORLD<T>& world_input) 
        :RENDERING_BLEND_SHADER<T>(shader1_input,shader2_input,direct_shading_only,world_input),field(field_input),low_value(low_value_input),high_value(high_value_input)
    {}

    T Blending_Fraction(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {const RENDERING_TRIANGULATED_SURFACE<T>& triangulated_surface=*(const RENDERING_TRIANGULATED_SURFACE<T>*)&intersection_object;
    assert(triangulated_surface.triangulated_surface.particles.array_collection->Size()==field.m);
    int i,j,k;triangulated_surface.triangulated_surface.mesh.elements(ray.ray.aggregate_id).Get(i,j,k);
    ARRAY_VIEW<const VECTOR<T,3> > X(triangulated_surface.triangulated_surface.particles.X);
    VECTOR<T,3> weights=TRIANGLE_3D<T>::Barycentric_Coordinates(triangulated_surface.Object_Space_Point(intersection_point),X(i),X(j),X(k));
    T blend=(weights.x*field(i)+weights.y*field(j)+weights.z*field(k)-low_value)/(high_value-low_value);
    return clamp(blend,(T)0,(T)1);}

//#####################################################################
};
}
#endif
