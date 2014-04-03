//#####################################################################
// Copyright 2003-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_TRANSPARENT_SHADER__
#define __RENDERING_TRANSPARENT_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_TRANSPARENT_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;
    T direct_reflection_coefficient_low,direct_reflection_coefficient_high;
    const RENDERING_BLEND_SHADER<T>* blend_shader;
    bool fresnel,ray_direction;

    RENDERING_TRANSPARENT_SHADER(const T constant_direct_reflection_coefficient,const bool fresnel_input,const bool ray_direction_input,RENDER_WORLD<T>& world_input)
        :MATERIAL_SHADER<T>(world_input),direct_reflection_coefficient_low(constant_direct_reflection_coefficient),direct_reflection_coefficient_high(0),
        blend_shader(0),fresnel(fresnel_input),ray_direction(ray_direction_input)
    {}

    RENDERING_TRANSPARENT_SHADER(const T direct_reflection_coefficient_low_input,const T direct_reflection_coefficient_high_input,const RENDERING_BLEND_SHADER<T>& blend_shader_input,
        const bool fresnel_input,const bool ray_direction_input,RENDER_WORLD<T>& world_input)
        :MATERIAL_SHADER<T>(world_input),direct_reflection_coefficient_low(direct_reflection_coefficient_low_input),direct_reflection_coefficient_high(direct_reflection_coefficient_high_input),
        blend_shader(&blend_shader_input),fresnel(fresnel_input),ray_direction(ray_direction_input)
    {}

    T Direct_Reflection_Coefficient(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {if(blend_shader){
        T fraction=blend_shader->Blending_Fraction(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
        return (1-fraction)*direct_reflection_coefficient_low+fraction*direct_reflection_coefficient_high;}
    else return direct_reflection_coefficient_low;}

    bool Create_Refraction_Ray(RENDERING_RAY<T>& transmitted_ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_RAY<T>& ray,const TV& same_side_normal,
        const TV& other_side_position,const TV& intersection_point) const
    {transmitted_ray=RENDERING_RAY<T>(RAY<TV>(other_side_position,ray.ray.direction,true),1,&entering_object);
    T exiting_index_of_refraction=exiting_object.Index_Of_Refraction(intersection_point);
    T entering_index_of_refraction=entering_object.Index_Of_Refraction(other_side_position);
    if(exiting_index_of_refraction!=entering_index_of_refraction){
        T ratio=exiting_index_of_refraction/entering_index_of_refraction;
        T determinant=sqr(ratio)*(sqr(TV::Dot_Product(same_side_normal,ray.ray.direction))-1)+1;
        if(determinant < 0)return true; // no transmitted ray
        else transmitted_ray.ray.direction=ratio*(ray.ray.direction-TV::Dot_Product(same_side_normal,ray.ray.direction)*same_side_normal)-sqrt(determinant)*same_side_normal;}
    return false;}

    void Compute_Coefficients(T& reflection_coefficient,T& transmission_coefficient,const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {T direct_reflection_coefficient=Direct_Reflection_Coefficient(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    if(fresnel){
        T R_dot_N=TV::Dot_Product(-ray.ray.direction,same_side_normal);
        reflection_coefficient=direct_reflection_coefficient+(1-direct_reflection_coefficient)*pow(1-R_dot_N,5);transmission_coefficient=1-reflection_coefficient;}
    else{reflection_coefficient=direct_reflection_coefficient;transmission_coefficient=1-reflection_coefficient;}}

    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {T reflection_coefficient,transmission_coefficient;
    Compute_Coefficients(reflection_coefficient,transmission_coefficient,ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    TV final_color;
    if(transmission_coefficient>0){
        TV other_side_position;
        if(ray_direction){other_side_position=intersection_point+ray.ray.direction*intersection_object.small_number*4;}
        else {other_side_position=intersection_point-same_side_normal*intersection_object.small_number*4;}
        RENDERING_RAY<T> refracted_ray;
        if(Create_Refraction_Ray(refracted_ray,exiting_object,entering_object,ray,same_side_normal,other_side_position,intersection_point)) reflection_coefficient=1;
        else{refracted_ray.ray_contribution=ray.ray_contribution*transmission_coefficient;final_color+=transmission_coefficient*world.Cast_Ray(refracted_ray,ray);}}
    if(reflection_coefficient>0){
        TV same_side_position;
        if(ray_direction){same_side_position=intersection_point-ray.ray.direction*intersection_object.small_number*2;}
        else{same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;}
        RENDERING_RAY<T> reflected_ray(RAY<TV>(same_side_position,world.Reflected_Direction(same_side_normal,ray.ray.direction),true),ray.ray_contribution*reflection_coefficient,&exiting_object);
        final_color+=reflection_coefficient*world.Cast_Ray(reflected_ray,ray);}
    return final_color;}

    TV Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal,const RENDERING_LIGHT<T>& light,const RENDERING_RAY<T>& full_ray) const
    {if(world.global_photon_map.photons.m || world.caustic_photon_map.photons.m) return TV();
    T reflection_coefficient,transmission_coefficient;
    Compute_Coefficients(reflection_coefficient,transmission_coefficient,ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    RENDERING_RAY<T> transmitted_ray;
    if(ray_direction){transmitted_ray=RENDERING_RAY<T>(RAY<TV>(SEGMENT_3D<T>(intersection_point+ray.ray.direction*intersection_object.small_number*10,full_ray.ray.Point(full_ray.ray.t_max))),                             ray.ray_contribution*transmission_coefficient,&entering_object);}
    else{transmitted_ray=RENDERING_RAY<T>(RAY<TV>(SEGMENT_3D<T>(intersection_point-same_side_normal*intersection_object.small_number*10,full_ray.ray.Point(full_ray.ray.t_max))),                             ray.ray_contribution*transmission_coefficient,&entering_object);}
    return transmission_coefficient*world.Incident_Light(transmitted_ray,light,full_ray,ray);}

    TV Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {T reflection_coefficient,transmission_coefficient;
    Compute_Coefficients(reflection_coefficient,transmission_coefficient,ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    TV final_color;
    if(transmission_coefficient>0){
        TV other_side_position=intersection_point-same_side_normal*intersection_object.small_number*2;
        RENDERING_RAY<T> refracted_ray;
        if(Create_Refraction_Ray(refracted_ray,exiting_object,entering_object,ray,same_side_normal,other_side_position,intersection_point)) reflection_coefficient=1;
        else{refracted_ray.ray_contribution=ray.ray_contribution*transmission_coefficient;final_color+=transmission_coefficient*world.Cast_Ray_For_Photon_Gather(refracted_ray,ray);}}
    if(reflection_coefficient>0){
        TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
        RENDERING_RAY<T> reflected_ray(RAY<TV>(same_side_position,world.Reflected_Direction(same_side_normal,ray.ray.direction),true),ray.ray_contribution*reflection_coefficient,&entering_object);
        final_color+=reflection_coefficient*world.Cast_Ray_For_Photon_Gather(reflected_ray,ray);}
        return final_color;}

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
                         const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal,const TV& power,
                         const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {T reflection_coefficient,transmission_coefficient;
    Compute_Coefficients(reflection_coefficient,transmission_coefficient,ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    T xi=world.random.Get_Uniform_Number((T)0,(T)1);
    if(xi<transmission_coefficient){
        // try to do transmission ray
        TV other_side_position=intersection_point-same_side_normal*intersection_object.small_number*4;
        RENDERING_RAY<T> refracted_ray;
        if(!Create_Refraction_Ray(refracted_ray,exiting_object,entering_object,ray,same_side_normal,other_side_position,intersection_point)){
            world.Cast_Photon(refracted_ray,ray,power,type,diffuse_bounces,specular_bounces+1);
            return;}}
    // reflect
    TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
    RENDERING_RAY<T> reflected_ray(RAY<TV>(same_side_position,world.Reflected_Direction(same_side_normal,ray.ray.direction),true),ray.ray_contribution*reflection_coefficient,&exiting_object);
    world.Cast_Photon(reflected_ray,ray,power,type,diffuse_bounces,specular_bounces+1);}

//#####################################################################
};
}
#endif
