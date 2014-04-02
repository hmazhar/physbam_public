//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_LAMBERTIAN_SHADER__
#define __RENDERING_LAMBERTIAN_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_LAMBERTIAN_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;

    const MATERIAL_SHADER<T>& shader;
    T ambient_coefficient,diffuse_coefficient;
    TV ambient_color,diffuse_color;
    TV ambient_factor,diffuse_factor; // precomputed multiplications
    bool visualize_photon_map_directly;

    RENDERING_LAMBERTIAN_SHADER(const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input,bool visualize_photon_map_directly_input=false) 
        :MATERIAL_SHADER<T>(world_input),shader(shader_input),visualize_photon_map_directly(visualize_photon_map_directly_input)
    {
        ambient_coefficient=0;diffuse_coefficient=1;ambient_color=TV(0,0,0);diffuse_color=TV(1,1,1);
        ambient_factor=ambient_color*ambient_coefficient;diffuse_factor=diffuse_color*diffuse_coefficient;
    }

    RENDERING_LAMBERTIAN_SHADER(const MATERIAL_SHADER<T>& shader_input,const T ambient_coefficient_input,const TV& ambient_color_input,
        const T diffuse_coefficient_input,const TV& diffuse_color_input,RENDER_WORLD<T>& world_input,bool visualize_photon_map_directly_input=false) 
        :MATERIAL_SHADER<T>(world_input),shader(shader_input),ambient_coefficient(ambient_coefficient_input),diffuse_coefficient(diffuse_coefficient_input),
        ambient_color(ambient_color_input),diffuse_color(diffuse_color_input),visualize_photon_map_directly(visualize_photon_map_directly_input)
    {
        ambient_factor=ambient_color*ambient_coefficient;diffuse_factor=diffuse_color*diffuse_coefficient;
    }

    TV Evaluate_Diffuse_BRDF(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {return diffuse_factor*shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal)/(T)pi;}

    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {if(visualize_photon_map_directly)return Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    const ARRAY<RENDERING_LIGHT<T> *>& lights=world.Lights();
    TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
    TV accumulated_color(0,0,0);
    for(int light_index=1;light_index<=lights.m;light_index++){
        if(lights(light_index)->photon_source_only) continue;
        TV accumulated_samples(0,0,0);
        ARRAY<RAY<TV> > sample_array;
        lights(light_index)->Sample_Points(same_side_position,same_side_normal,sample_array);
        for(int sample=1;sample<=sample_array.m;sample++){
            RENDERING_RAY<T> ray_to_light(sample_array(sample),1,&exiting_object);
            TV light_color=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,ray);
            T L_N=TV::Dot_Product(ray_to_light.ray.direction,same_side_normal);
            if(L_N<0) continue;
            accumulated_samples+=L_N*light_color;}
        accumulated_color+=accumulated_samples/T(sample_array.m);}
    return (accumulated_color+ambient_factor)*Evaluate_Diffuse_BRDF(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    TV Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {if(visualize_photon_map_directly)return TV(0,0,0);// don't want to double count caustics
    TV caustic_light=world.caustic_photon_map.Irradiance_Estimate(intersection_point,same_side_normal,world.max_photon_distance,
        world.number_of_photons_for_estimate,ray,PHOTON_MAP<T>::CAUSTIC_PHOTON_MAP);
    TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
    TV indirect_irradiance;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(world.global_photon_map.photons.m){//indirect_irradiance=world.irradiance_cache.Compute_Indirect_Light(world,ray,exiting_object,entering_object,intersection_object,same_side_position,same_side_normal);
        if(world.use_irradiance_cache)indirect_irradiance=world.irradiance_cache.Compute_Indirect_Light(world,ray,exiting_object,entering_object,intersection_object,same_side_position,same_side_normal);
        else indirect_irradiance=world.global_photon_map.Irradiance_Estimate(same_side_position,same_side_normal,world.max_photon_distance,
            world.number_of_photons_for_estimate,ray,PHOTON_MAP<T>::GLOBAL_PHOTON_MAP);}
#endif
    return (caustic_light+indirect_irradiance)*Evaluate_Diffuse_BRDF(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    TV Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {TV approximate_color=world.global_photon_map.Irradiance_Estimate(intersection_point,same_side_normal,world.max_photon_distance,
            world.number_of_photons_for_estimate,ray,PHOTON_MAP<T>::GLOBAL_PHOTON_MAP);
    if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Approx color %f %f %f\n",approximate_color.x,approximate_color.y,approximate_color.z));
    return approximate_color*Evaluate_Diffuse_BRDF(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal,const TV& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {TV reflectance=diffuse_factor*shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    T maximum_power=power.Max();
    TV reflectance_times_power=reflectance*power;
    T probability_diffuse=reflectance_times_power.Max()/maximum_power;
    T xi=world.random.Get_Uniform_Number(T(0),T(1));
    if(xi<probability_diffuse){
        if(type==PHOTON_MAP<T>::GLOBAL_PHOTON_MAP){ //only do this if we are photon tracing for global map (only way we'd store diffuse-diffuse photons)
            TV reflected_photon_power=reflectance_times_power/probability_diffuse;
            TV reflected_direction=world.Random_Reflected_Direction_In_Hemisphere(same_side_normal,world.random.Get_Uniform_Number(T(0),T(1)),world.random.Get_Uniform_Number(T(0),T(1)));
            TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
            RENDERING_RAY<T> photon_ray(RAY<TV>(same_side_position,reflected_direction),ray.ray_contribution,&exiting_object);
            world.Cast_Photon(photon_ray,ray,reflected_photon_power,type,diffuse_bounces+1,specular_bounces);}}
    else{
        if(type==PHOTON_MAP<T>::CAUSTIC_PHOTON_MAP&&specular_bounces>0&&diffuse_bounces==0)
            world.caustic_photon_map.Store_Photon(intersection_point,ray.ray.direction,power,ray.recursion_depth);
        else if(type==PHOTON_MAP<T>::GLOBAL_PHOTON_MAP) 
            world.global_photon_map.Store_Photon(intersection_point,ray.ray.direction,power,ray.recursion_depth);}}
    
//#####################################################################
};
}
#endif
