//#####################################################################
// Copyright 2003-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_COLOR_BLEND_SHADER__
#define __RENDERING_COLOR_BLEND_SHADER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_COLOR_BLEND_SHADER:public MATERIAL_SHADER<T>
{
public:
    using MATERIAL_SHADER<T>::world;

    ARRAY<const MATERIAL_SHADER<T>*> shaders;
    ARRAY<T> weights;
    T total_weight;

    RENDERING_COLOR_BLEND_SHADER(RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),total_weight(0)
    {}

    void Add_Shader(const MATERIAL_SHADER<T>* shader,T weight)
    {shaders.Append(shader);weights.Append(weight);total_weight+=weight;}

    VECTOR<T,3> Evaluate_Diffuse_BRDF(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> color_accumulator;
    for(int i=1;i<=shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Evaluate_Diffuse_BRDF(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return color_accumulator;}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> color_accumulator;
    for(int i=1;i<=shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal); 
    return color_accumulator;}

    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,
        const RENDERING_RAY<T>& full_ray) const
    {VECTOR<T,3> color_accumulator;
    for(int i=1;i<=shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Light_Ray(ray,exiting_object,entering_object,intersection_object,intersection_point,
        same_side_normal,light,full_ray);
    return color_accumulator;}

    VECTOR<T,3> Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> color_accumulator;
    for(int i=1;i<=shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Surface_Using_Indirect_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal); 
    return color_accumulator;}

    VECTOR<T,3> Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> color_accumulator;
    for(int i=1;i<=shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return color_accumulator;}

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {T r=world.random.Get_Uniform_Number((T)0,total_weight);
    for(int i=1;i<=shaders.m-1;i++){
        if(r<=weights(i)){shaders(i)->Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);return;}
        else r-=weights(i);}
    shaders(shaders.m)->Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);}

//#####################################################################
};
}
#endif
