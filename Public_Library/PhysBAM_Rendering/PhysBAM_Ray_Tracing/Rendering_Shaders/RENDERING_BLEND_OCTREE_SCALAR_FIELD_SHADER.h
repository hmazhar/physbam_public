#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_BLEND_OCTREE_SCALAR_FIELD_SHADER__
#define __RENDERING_BLEND_OCTREE_SCALAR_FIELD_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OCTREE_IMPLICIT_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

// mostly taken from RENDERING_TRANSPARENT_SHADER
template<class T>
class RENDERING_BLEND_OCTREE_SCALAR_FIELD_SHADER:public MATERIAL_SHADER<T>
{
public:
    using MATERIAL_SHADER<T>::world;
    T direct_reflection_coefficient1,direct_reflection_coefficient2;
    bool fresnel, ray_direction;

    const ARRAY<T>* values1;
    const ARRAY<T>* values2; // may be null
    bool use_abs_value;
    T tolerance;
    T blend_band;
    int blend_mode;

    RENDERING_BLEND_OCTREE_SCALAR_FIELD_SHADER(T direct_reflection_coefficient1_input,T direct_reflection_coefficient2_input,
                                               bool fresnel_input,const ARRAY<T>* values1_input,const ARRAY<T>* values2_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),direct_reflection_coefficient1(direct_reflection_coefficient1_input),direct_reflection_coefficient2(direct_reflection_coefficient2_input),
        fresnel(fresnel_input),ray_direction(false),values1(values1_input),values2(values2_input),use_abs_value(false),tolerance((T)1e-8),blend_band((T).1),blend_mode(1)
    {}

    bool Create_Refraction_Ray(RENDERING_RAY<T>& transmitted_ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_RAY<T>& ray,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& other_side_position) const
    {transmitted_ray=RENDERING_RAY<T>(RAY<VECTOR<T,3> >(other_side_position,ray.ray.direction,true),1,&entering_object);
    if(exiting_object.index_of_refraction!=entering_object.index_of_refraction){
        T ratio=exiting_object.index_of_refraction/entering_object.index_of_refraction;
        T determinant=sqr(ratio)*(sqr(VECTOR<T,3>::Dot_Product(same_side_normal,ray.ray.direction))-1)+1;
        if(determinant < 0)return true; // no transmitted ray
        else transmitted_ray.ray.direction=ratio*(ray.ray.direction-VECTOR<T,3>::Dot_Product(same_side_normal,ray.ray.direction)*same_side_normal)-sqrt(determinant)*same_side_normal;}
    return false;}

    void Compute_Coefficients(const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,T& reflection_coefficient,T& transmission_coefficient,const RENDERING_RAY<T>& ray,const VECTOR<T,3>& same_side_normal) const
    {const OCTREE_GRID<T>& grid=((const RENDERING_OCTREE_IMPLICIT_SURFACE<T>*)&intersection_object)->implicit_surface->levelset.grid;
    const OCTREE_CELL<T>* cell=grid.Leaf_Cell(intersection_point);
    T value1=LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,cell,*values1,intersection_point);
    T blend=0;
    if(blend_mode==1){assert(values2);
        T value2=LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,cell,*values2,intersection_point);
        if(use_abs_value){value1=abs(value1);value2=abs(value2);}
        blend=clamp((T)0.5+(value1-value2)/blend_band,(T)0,(T)1);}
    else{assert(blend_mode==2);
        blend=clamp(value1/blend_band,(T)0,(T)1);}

    T direct_reflection_coefficient=(1-blend)*direct_reflection_coefficient1+blend*direct_reflection_coefficient2;
    if(fresnel){
        T R_dot_N=VECTOR<T,3>::Dot_Product(-ray.ray.direction,same_side_normal);
        reflection_coefficient=direct_reflection_coefficient+(1-direct_reflection_coefficient)*pow(1-R_dot_N,5);transmission_coefficient=1-reflection_coefficient;}
    else{reflection_coefficient=direct_reflection_coefficient;transmission_coefficient=1-reflection_coefficient;}}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {T reflection_coefficient,transmission_coefficient;Compute_Coefficients(intersection_object,intersection_point,reflection_coefficient,transmission_coefficient,ray,same_side_normal);
    VECTOR<T,3> final_color;
    if(transmission_coefficient>0){
        VECTOR<T,3> other_side_position;
        if (ray_direction){other_side_position=intersection_point+ray.ray.direction*intersection_object.small_number*4;}
        else {other_side_position=intersection_point-same_side_normal*intersection_object.small_number*4;}
        RENDERING_RAY<T> refracted_ray; 
        if(Create_Refraction_Ray(refracted_ray,exiting_object,entering_object,ray,same_side_normal,other_side_position)) reflection_coefficient=1;
        else{refracted_ray.ray_contribution=ray.ray_contribution*transmission_coefficient;final_color+=transmission_coefficient*world.Cast_Ray(refracted_ray,ray);}}
    if(reflection_coefficient>0){
        VECTOR<T,3> same_side_position;
        if(ray_direction){same_side_position=intersection_point-ray.ray.direction*intersection_object.small_number*2;}
        else{same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;}
        RENDERING_RAY<T> reflected_ray(RAY<VECTOR<T,3> >(same_side_position,world.Reflected_Direction(same_side_normal,ray.ray.direction),true),ray.ray_contribution*reflection_coefficient,&exiting_object);
        final_color+=reflection_coefficient*world.Cast_Ray(reflected_ray,ray);}
    return final_color;}

    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,const RENDERING_RAY<T>& full_ray) const
    {if(world.global_photon_map.photons.m)return VECTOR<T,3>();
    T reflection_coefficient,transmission_coefficient;Compute_Coefficients(intersection_object,intersection_point,reflection_coefficient,transmission_coefficient,ray,same_side_normal);
    RENDERING_RAY<T> transmitted_ray;
    if (ray_direction){transmitted_ray=RENDERING_RAY<T>(RAY<VECTOR<T,3> >(SEGMENT_3D<T>(intersection_point+ray.ray.direction*intersection_object.small_number*10,full_ray.ray.Point(full_ray.ray.t_max))),                             ray.ray_contribution*transmission_coefficient,&entering_object);}
    else{transmitted_ray=RENDERING_RAY<T>(RAY<VECTOR<T,3> >(SEGMENT_3D<T>(intersection_point-same_side_normal*intersection_object.small_number*10,full_ray.ray.Point(full_ray.ray.t_max))),                             ray.ray_contribution*transmission_coefficient,&entering_object);}
    return transmission_coefficient*world.Incident_Light(transmitted_ray,light,full_ray,ray);}

    VECTOR<T,3> Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {T reflection_coefficient,transmission_coefficient;Compute_Coefficients(intersection_object,intersection_point,reflection_coefficient,transmission_coefficient,ray,same_side_normal);
    VECTOR<T,3> final_color;
    if(transmission_coefficient>0){
        VECTOR<T,3> other_side_position=intersection_point-same_side_normal*intersection_object.small_number*2;
        RENDERING_RAY<T> refracted_ray; 
        if(Create_Refraction_Ray(refracted_ray,exiting_object,entering_object,ray,same_side_normal,other_side_position)) reflection_coefficient=1;
        else{refracted_ray.ray_contribution=ray.ray_contribution*transmission_coefficient;final_color+=transmission_coefficient*world.Cast_Ray_For_Photon_Gather(refracted_ray,ray);}}
    if(reflection_coefficient>0){
        VECTOR<T,3> same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
        RENDERING_RAY<T> reflected_ray(RAY<VECTOR<T,3> >(same_side_position,world.Reflected_Direction(same_side_normal,ray.ray.direction),true),ray.ray_contribution*reflection_coefficient,&entering_object);
        final_color+=reflection_coefficient*world.Cast_Ray_For_Photon_Gather(reflected_ray,ray);}
        return final_color;}

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
                         const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
                         const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
};
}
#endif
#endif
