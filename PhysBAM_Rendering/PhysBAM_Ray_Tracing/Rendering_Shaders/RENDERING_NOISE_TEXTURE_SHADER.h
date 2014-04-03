//#####################################################################
// Copyright 2004-2007, Igor Neverov, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_NOISE_TEXTURE_SHADER  
//#####################################################################
#ifndef __RENDERING_NOISE_TEXTURE_SHADER__
#define __RENDERING_NOISE_TEXTURE_SHADER__

#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_NOISE_TEXTURE_SHADER:public MATERIAL_SHADER<T>
{
public:
    VECTOR<T,3> color0,color1;
    T texture_scaling_factor;

    RENDERING_NOISE_TEXTURE_SHADER(const VECTOR<T,3>& color0_input,const VECTOR<T,3>& color1_input,const T texture_scaling_factor_input,RENDER_WORLD<T>& world) 
        :MATERIAL_SHADER<T>(world),color0(color0_input),color1(color1_input),texture_scaling_factor(texture_scaling_factor_input)
    {}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {RAY<VECTOR<T,3> > object_space_ray=intersection_object.Object_Space_Ray(ray.ray);VECTOR<T,3> position=object_space_ray.Point(object_space_ray.t_max);
     T s,t;intersection_object.Get_Texture_Coordinates(position,ray.ray.aggregate_id,s,t);VECTOR<double,3> p(s*texture_scaling_factor,t*texture_scaling_factor,0);
     T a=(T)NOISE<double>::Noise1(p,10,(T).5)/4;a=clamp<T>((a-(T).17)*(T)1.5625,0,1);return (1-a)*color0+a*color1;}

//#####################################################################
};
}
#endif

