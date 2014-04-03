//#####################################################################
// Copyright 2004, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_PHONG_SHADER__
#define __RENDERING_PHONG_SHADER__

#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_PHONG_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;

    TV ambient_color,diffuse_color,specular_color;
    T specular_exponent;
    const MATERIAL_SHADER<T>& shader;

    RENDERING_PHONG_SHADER(const TV& ambient_color_input,const TV& diffuse_color_input,const TV& specular_color_input,const T specular_exponent_input,
        const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),ambient_color(ambient_color_input),diffuse_color(diffuse_color_input),specular_color(specular_color_input),specular_exponent(specular_exponent_input),
         shader(shader_input)
    {}
    
    RENDERING_PHONG_SHADER(const TV& ambient_color_input,const TV& diffuse_color_input,const T diffuse_coefficient,const TV& specular_color_input,
        const T specular_coefficient,const T specular_exponent_input,const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),ambient_color(ambient_color_input),diffuse_color(diffuse_color_input*diffuse_coefficient),specular_color(specular_color_input*specular_coefficient),
        specular_exponent(specular_exponent_input),shader(shader_input)
    {}
    
    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {const ARRAY<RENDERING_LIGHT<T>*>& lights=world.Lights();
    TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number,reflected_direction=ray.ray.Reflected_Direction(same_side_normal);
    TV accumulated_diffuse_color(0,0,0),accumulated_specular_color(0,0,0);
    for(int light_index=1;light_index<=lights.m;light_index++){
        TV accumulated_diffuse_samples(0,0,0),accumulated_specular_samples(0,0,0);
        ARRAY<RAY<TV> > sample_array;
        lights(light_index)->Sample_Points(same_side_position,same_side_normal,sample_array);
        for(int sample=1;sample<=sample_array.m;sample++){
            RENDERING_RAY<T> ray_to_light(sample_array(sample),1,&entering_object);
            TV light_color=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,ray);
            T L_N=TV::Dot_Product(ray_to_light.ray.direction,same_side_normal);
            if(L_N<0) continue;
            accumulated_diffuse_samples+=L_N*light_color;
            accumulated_specular_samples+=std::pow(PhysBAM::max(TV::Dot_Product(ray_to_light.ray.direction,reflected_direction),T(0)),specular_exponent)*light_color;}
        accumulated_diffuse_color+=accumulated_diffuse_samples/T(sample_array.m);
        accumulated_specular_color+=accumulated_specular_samples/T(sample_array.m);}
    TV ret=(ambient_color+accumulated_diffuse_color*diffuse_color+accumulated_specular_color*specular_color)*
        shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return ret;}

//#####################################################################
};
}
#endif
