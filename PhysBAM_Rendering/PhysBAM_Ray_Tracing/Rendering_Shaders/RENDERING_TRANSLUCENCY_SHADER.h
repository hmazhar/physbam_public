//#####################################################################
// Copyright 2004, Frank Losasso, Andrew Selle, Jiayi Chong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_TRANSLUCENCY_SHADER__
#define __RENDERING_TRANSLUCENCY_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_TRANSLUCENCY_SHADER:public MATERIAL_SHADER<T>
{
public:

    using MATERIAL_SHADER<T>::world;

    const MATERIAL_SHADER<T>& shader;
    T translucency;

    RENDERING_TRANSLUCENCY_SHADER(const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input,T translucency_input) 
        :MATERIAL_SHADER<T>(world_input),shader(shader_input),translucency(translucency_input)
    {}

    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,const RENDERING_RAY<T>& full_ray) const
    {if(world.global_photon_map.photons.m)return VECTOR<T,3>();
    RENDERING_RAY<T> transmitted_ray=RENDERING_RAY<T>(RAY<VECTOR<T,3> >(SEGMENT_3D<T>(intersection_point+ray.ray.direction*intersection_object.small_number*3,full_ray.ray.Point(full_ray.ray.t_max))),
        ray.ray_contribution,&entering_object);
    return world.Incident_Light(transmitted_ray,light,full_ray,ray);}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> other_side_position=intersection_point-same_side_normal*intersection_object.small_number*4;
    RENDERING_RAY<T> transmitted_ray(RAY<VECTOR<T,3> >(other_side_position,ray.ray.direction,true),1,&entering_object);
    return (1-translucency)*world.Cast_Ray(transmitted_ray,ray)+translucency*shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,
        intersection_point,same_side_normal);}

//#####################################################################
};
}
#endif
