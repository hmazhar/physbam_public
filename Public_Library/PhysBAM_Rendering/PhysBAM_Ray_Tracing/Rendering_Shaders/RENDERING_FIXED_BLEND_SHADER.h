//#####################################################################
// Copyright 2008, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_FIXED_BLEND_SHADER__
#define __RENDERING_FIXED_BLEND_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_FIXED_BLEND_SHADER:public RENDERING_BLEND_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    T blending_fraction;
    
    RENDERING_FIXED_BLEND_SHADER(const T blending_fraction_input,const MATERIAL_SHADER<T>& shader1_input,const MATERIAL_SHADER<T>& shader2_input,const bool direct_shading_only,
        RENDER_WORLD<T>& world_input)
        :RENDERING_BLEND_SHADER<T>(shader1_input,shader2_input,direct_shading_only,world_input),blending_fraction(blending_fraction_input)
    {}

    T Blending_Fraction(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal) const
    {return blending_fraction;}

//#####################################################################
};
}
#endif
