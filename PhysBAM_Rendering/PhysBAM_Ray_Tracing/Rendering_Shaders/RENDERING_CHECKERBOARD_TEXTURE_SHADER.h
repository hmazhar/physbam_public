//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_CHECKERBOARD_TEXTURE_SHADER__
#define __RENDERING_CHECKERBOARD_TEXTURE_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_CHECKERBOARD_TEXTURE_SHADER:public MATERIAL_SHADER<T>
{
public:
    const MATERIAL_SHADER<T> &shader1,&shader2;
    T texture_scaling_factor;

    RENDERING_CHECKERBOARD_TEXTURE_SHADER(const MATERIAL_SHADER<T>& shader1_input,const MATERIAL_SHADER<T>& shader2_input,const T texture_scaling_factor_input,
        RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),shader1(shader1_input),shader2(shader2_input),texture_scaling_factor(texture_scaling_factor_input)
    {}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> position=ray.ray.Point(ray.ray.t_max);T s,t;intersection_object.Get_Texture_Coordinates(position,ray.ray.aggregate_id,s,t);
    if((((int)floor(s*texture_scaling_factor)&0x1)==0)^(((int)floor(t*texture_scaling_factor)&0x1)==0))
        return shader1.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    else return shader2.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);}

//#####################################################################
};
}
#endif
