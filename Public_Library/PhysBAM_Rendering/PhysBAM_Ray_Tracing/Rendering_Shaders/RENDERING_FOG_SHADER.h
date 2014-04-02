//#####################################################################
// Copyright 2003, 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_FOG_SHADER__
#define __RENDERING_FOG_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_FOG_SHADER:public VOLUMETRIC_SHADER<T>
{
public:
    T density;
    VECTOR<T,3> fog_color;

    RENDERING_FOG_SHADER(const T density_input,const VECTOR<T,3>& fog_color_input,RENDER_WORLD<T>& world_input) 
        :VOLUMETRIC_SHADER<T>(world_input),density(density_input),fog_color(fog_color_input)
    {}

    VECTOR<T,3> Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const VECTOR<T,3>& color) PHYSBAM_OVERRIDE
    {if(ray.ray.semi_infinite) return VECTOR<T,3>(fog_color);
    T alpha=clamp(exp(-density*ray.ray.t_max),(T)0,(T)1); 
    return color*alpha+fog_color*(1-alpha);}

    VECTOR<T,3> Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,const VECTOR<T,3>& light_color) PHYSBAM_OVERRIDE
    {return Attenuate_Color(ray,object,light_color);}

//#####################################################################
};
}
#endif
