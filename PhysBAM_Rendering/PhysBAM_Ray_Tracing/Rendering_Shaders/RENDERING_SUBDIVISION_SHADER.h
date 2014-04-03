//#####################################################################
// Copyright 2008, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_SUBDIVISION_SHADER__
#define __RENDERING_SUBDIVISION_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_SUBDIVISION_SHADER:public MATERIAL_SHADER<T>
{
public:

    using MATERIAL_SHADER<T>::world;

    const MATERIAL_SHADER<T>& shader;

    RENDERING_SUBDIVISION_SHADER(const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),shader(shader_input)
    {}

    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,const RENDERING_RAY<T>& full_ray) const
    {return shader.Shade_Light_Ray(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,light,full_ray);}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {// recast ray from start, using secondary intersection for object and requiring a hit.
        RENDERING_RAY<T> copy_ray(ray);
        copy_ray.ray.t_max=(T)0;copy_ray.ray.semi_infinite=true;
        const RENDERING_TRIANGULATED_SURFACE<T>* surface=dynamic_cast<const RENDERING_TRIANGULATED_SURFACE<T>*>(&intersection_object);
        if(!surface) return VECTOR<T,3>();
        TRIANGULATED_SURFACE<T>* temp=surface->base_triangulated_surface;
        surface->base_triangulated_surface=&surface->triangulated_surface;
        if(surface->Intersection(copy_ray.ray)){
            VECTOR<T,3> intersection_point=copy_ray.ray.Point(copy_ray.ray.t_max);
            bool flipped_normal;VECTOR<T,3> same_side_normal=world.Same_Side_Normal(copy_ray,intersection_point,intersection_object,flipped_normal);
            VECTOR<T,3> color=shader.Shade_Surface_Using_Direct_Illumination(copy_ray,*ray.current_object,entering_object,intersection_object,intersection_point,same_side_normal);
            surface->base_triangulated_surface=temp;
            return color;}

        surface->base_triangulated_surface=temp;
        // suppress this object and recast
        int priority=exiting_object.priority;
        intersection_object.priority=-2;
        copy_ray.ray.t_max=(T)0;copy_ray.ray.semi_infinite=true;
        VECTOR<T,3> color=world.Cast_Ray(copy_ray,ray);
        intersection_object.priority=priority;
        return color;}
        
//#####################################################################
};
}
#endif
