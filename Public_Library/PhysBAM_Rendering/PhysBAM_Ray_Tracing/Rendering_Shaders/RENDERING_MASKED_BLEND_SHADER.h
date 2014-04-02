//#####################################################################
// Copyright 2008, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_MASKED_BLEND_SHADER__
#define __RENDERING_MASKED_BLEND_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_MASKED_BLEND_SHADER:public RENDERING_BLEND_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    const MATERIAL_SHADER<T> &mask_shader;
    int mask_usage; // 0=average, 1=red channel, 2=green channel, 3=blue channel

    RENDERING_MASKED_BLEND_SHADER(const MATERIAL_SHADER<T>& mask_shader_input,const MATERIAL_SHADER<T>& shader1_input,const MATERIAL_SHADER<T>& shader2_input,
        const bool direct_shading_only_input,RENDER_WORLD<T>& world_input)
        :RENDERING_BLEND_SHADER<T>(shader1_input,shader2_input,direct_shading_only_input,world_input),mask_shader(mask_shader_input),mask_usage(3)
    {}

    void Use_Gray()
    {mask_usage=0;}

    void Use_Red()
    {mask_usage=1;}

    void Use_Green()
    {mask_usage=2;}

    void Use_Blue()
    {mask_usage=3;}

    T Blending_Fraction(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const PHYSBAM_OVERRIDE
    {TV mask=mask_shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return mask_usage?mask(mask_usage):mask.Average();}
//#####################################################################
};
}
#endif
