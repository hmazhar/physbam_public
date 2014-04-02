//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_SHELL_EMISSION_SHADER__
#define __RENDERING_SHELL_EMISSION_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_SHELL_EMISSION_SHADER:public MATERIAL_SHADER<T>
{
public:
    using MATERIAL_SHADER<T>::world;

    T shell_thickness;
    T shell_radius_of_curvature;
    T shell_amplification_factor;
    VECTOR<T,3> shell_emission_color;
    
    RENDERING_SHELL_EMISSION_SHADER(const T shell_thickness,const T shell_radius_of_curvature,const T shell_amplification_factor,const VECTOR<T,3>& shell_emission_color,
        RENDER_WORLD<T>& world_input)
        :MATERIAL_SHADER<T>(world_input),shell_thickness(shell_thickness),shell_radius_of_curvature(shell_radius_of_curvature),shell_amplification_factor(shell_amplification_factor),
        shell_emission_color(shell_emission_color)
    {}

    VECTOR<T,3> Shell_Emission(const VECTOR<T,3>& ray_direction,const VECTOR<T,3>& normal)const
    {T y=shell_radius_of_curvature*abs(VECTOR<T,3>::Dot_Product(ray_direction,normal)),x=sqrt(sqr(shell_radius_of_curvature)-sqr(y));
    T inner_radius=shell_radius_of_curvature-shell_thickness;
    if(x>inner_radius) return y*shell_emission_color*shell_amplification_factor; // doesn't intersect the inner shell
    else{T y_inside=sqrt(sqr(inner_radius)-sqr(x));return (y-y_inside)*shell_emission_color*shell_amplification_factor;}}

    void Create_Transmission_Ray(RENDERING_RAY<T>& transmitted_ray,const RENDERING_OBJECT<T>& entering_object,const RENDERING_RAY<T>& ray,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> other_side_position=intersection_point+ray.ray.direction*intersection_object.small_number*4;
    transmitted_ray=RENDERING_RAY<T>(RAY<VECTOR<T,3> >(other_side_position,ray.ray.direction,true),1,&entering_object);}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {RENDERING_RAY<T> transmitted_ray;
    Create_Transmission_Ray(transmitted_ray,entering_object,ray,intersection_object,intersection_point,same_side_normal);
    return world.Cast_Ray(transmitted_ray,ray)+Shell_Emission(ray.ray.direction,same_side_normal);}

    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,
        const RENDERING_RAY<T>& full_ray) const
    {if(world.use_photon_mapping)return VECTOR<T,3>(0,0,0);
    RENDERING_RAY<T> transmitted_ray;
    Create_Transmission_Ray(transmitted_ray,entering_object,ray,intersection_object,intersection_point,same_side_normal);
    return world.Incident_Light(transmitted_ray,light,full_ray,ray);}

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {RENDERING_RAY<T> transmitted_ray; 
    Create_Transmission_Ray(transmitted_ray,entering_object,ray,intersection_object,intersection_point,same_side_normal);
    world.Cast_Photon(transmitted_ray,ray,power,type,diffuse_bounces,specular_bounces+1);}

//#####################################################################
};
}
#endif
