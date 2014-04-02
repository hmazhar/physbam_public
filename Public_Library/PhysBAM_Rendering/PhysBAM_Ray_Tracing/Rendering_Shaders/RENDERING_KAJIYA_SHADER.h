//#####################################################################
// Copyright 2007, Joyce Pan.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_KAJIYA_SHADER__
#define __RENDERING_KAJIYA_SHADER__

#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <cmath>
namespace PhysBAM{

template<class T>
class RENDERING_KAJIYA_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;

    TV diffuse_color,specular_color;
    T diffuse_coefficient,specular_coefficient,specular_exponent;
    
    RENDERING_KAJIYA_SHADER(const TV& diffuse_color_input,const T diffuse_coefficient_input,const TV& specular_color_input,
        const T specular_coefficient_input,const T specular_exponent_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),diffuse_color(diffuse_color_input),specular_color(specular_color_input),diffuse_coefficient(diffuse_coefficient_input),
        specular_coefficient(specular_coefficient_input),specular_exponent(specular_exponent_input)
    {}
    
    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {const ARRAY<RENDERING_LIGHT<T>*>& lights=world.Lights();
    TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number;
    TV world_tangent,world_bitangent;
    intersection_object.Get_World_Space_Tangent_And_Bitangent(intersection_point,same_side_normal,ray.ray.aggregate_id,world_tangent,world_bitangent);
    TV accumulated_diffuse_color,accumulated_specular_color;
    for(int light_index=1;light_index<=lights.m;light_index++){
        TV accumulated_diffuse_samples,accumulated_specular_samples;
        ARRAY<RAY<VECTOR<T,3> > > sample_array;
        lights(light_index)->Sample_Points(same_side_position,same_side_normal,sample_array);
        for(int sample=1;sample<=sample_array.m;sample++){
            RENDERING_RAY<T> ray_to_light(sample_array(sample),1,&entering_object);
            //TV light_color=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,ray);
            TV light_color=lights(light_index)->Emitted_Light(ray_to_light);
            T cos_T_L=TV::Dot_Product(world_tangent,ray_to_light.ray.direction),cos_T_V=TV::Dot_Product(world_tangent,-ray.ray.direction);
            T theta_T_L=acos(cos_T_L),theta_T_V=acos(cos_T_V);
            T sin_T_L=sin(theta_T_L),sin_T_V=sin(theta_T_V);
            accumulated_diffuse_samples+=sin_T_L*light_color;
            accumulated_specular_samples+=std::pow(cos_T_L*cos_T_V+sin_T_L*sin_T_V,specular_exponent)*light_color;}
        accumulated_diffuse_color+=accumulated_diffuse_samples/T(sample_array.m);
        accumulated_specular_color+=accumulated_specular_samples/T(sample_array.m);}
    return accumulated_diffuse_color*diffuse_color*diffuse_coefficient+accumulated_specular_color*specular_color*specular_coefficient;} 

//#####################################################################
};
}
#endif
