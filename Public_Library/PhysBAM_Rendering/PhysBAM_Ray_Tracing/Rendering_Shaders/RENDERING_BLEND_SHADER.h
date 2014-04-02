//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_BLEND_SHADER__
#define __RENDERING_BLEND_SHADER__

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BLEND_SHADER:public MATERIAL_SHADER<T>
{
public:
    using MATERIAL_SHADER<T>::world;

    const MATERIAL_SHADER<T> &shader1,&shader2;

    // TODO: remove siggraph efficiency hack
    const bool direct_shading_only;
    
    RENDERING_BLEND_SHADER(const MATERIAL_SHADER<T>& shader1_input,const MATERIAL_SHADER<T>& shader2_input,const bool direct_shading_only_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),shader1(shader1_input),shader2(shader2_input),direct_shading_only(direct_shading_only_input)
    {}

    VECTOR<T,3> Evaluate_Diffuse_BRDF(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {T weight=Blending_Fraction(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return (1-weight)*shader1.Evaluate_Diffuse_BRDF(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal)
        +weight*shader2.Evaluate_Diffuse_BRDF(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {T weight=Blending_Fraction(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return (1-weight)*shader1.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal)
        +weight*shader2.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,
        const RENDERING_RAY<T>& full_ray) const
    {if(direct_shading_only) return VECTOR<T,3>(0,0,0);
    T weight=Blending_Fraction(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return (1-weight)*shader1.Shade_Light_Ray(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,light,full_ray);
    +weight*shader2.Shade_Light_Ray(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,light,full_ray);}

    VECTOR<T,3> Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {if(direct_shading_only) return VECTOR<T,3>(0,0,0);
    T weight=Blending_Fraction(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return (1-weight)*shader1.Shade_Surface_Using_Indirect_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal)
        +weight*shader2.Shade_Surface_Using_Indirect_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    VECTOR<T,3> Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {if(direct_shading_only) return VECTOR<T,3>(0,0,0);
    T weight=Blending_Fraction(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return (1-weight)*shader1.Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal)
        +weight*shader2.Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}
    
    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {T weight=Blending_Fraction(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    if(world.random.Get_Uniform_Number((T)0,(T)1)>weight) shader1.Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);
    else shader2.Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);}

//#####################################################################
    virtual T Blending_Fraction(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const=0;
//#####################################################################
};
}
#endif
