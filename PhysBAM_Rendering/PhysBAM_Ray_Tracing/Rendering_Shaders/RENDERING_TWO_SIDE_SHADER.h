//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_TWO_SIDE_SHADER__
#define __RENDERING_TWO_SIDE_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_TWO_SIDE_SHADER:public MATERIAL_SHADER<T>
{
public:
    const MATERIAL_SHADER<T> &front_shader,&back_shader;

    RENDERING_TWO_SIDE_SHADER(const MATERIAL_SHADER<T> &front_shader_input,const MATERIAL_SHADER<T> &back_shader_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),front_shader(front_shader_input),back_shader(back_shader_input)
    {}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {return VECTOR<T,3>::Dot_Product(ray.ray.direction,intersection_object.Normal(intersection_point,ray.ray.aggregate_id))>0?
        back_shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal):
        front_shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    VECTOR<T,3> Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {return VECTOR<T,3>::Dot_Product(ray.ray.direction,intersection_object.Normal(intersection_point,ray.ray.aggregate_id))>0?
        back_shader.Shade_Surface_Using_Indirect_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal):
        front_shader.Shade_Surface_Using_Indirect_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

    VECTOR<T,3> Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {return VECTOR<T,3>::Dot_Product(ray.ray.direction,intersection_object.Normal(intersection_point,ray.ray.aggregate_id))>0?
        back_shader.Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal):
        front_shader.Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}
    
    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {if(VECTOR<T,3>::Dot_Product(ray.ray.direction,intersection_object.Normal(intersection_point,ray.ray.aggregate_id))>0)
        back_shader.Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);
    else front_shader.Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);}

//#####################################################################
};
}
#endif
