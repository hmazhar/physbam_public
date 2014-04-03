//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_MARBLE_SHADER  
//#####################################################################
#ifndef __RENDERING_MARBLE_SHADER__
#define __RENDERING_MARBLE_SHADER__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <iostream>
namespace PhysBAM{

template<class T>
class RENDERING_MARBLE_SHADER:public MATERIAL_SHADER<T>
{
public:
    typedef VECTOR<T,3> TV;

    VECTOR<T,3> color1,color2;
    int octaves;
    T lacunarity,gain;
    T low,high;
    T vain_value,vain_width,clamp_width; // vain is dark up to width and then further width from the vain clamps to fullest extent
    INTERPOLATION_CURVE<T,T> curve;

    T Smooth(const T value,const T min,const T max) const
    {T v=clamp((value-min)/(max-min),(T)0,(T)1);
    return v*v*(3-2*value);}

    RENDERING_MARBLE_SHADER(const VECTOR<T,3>& color1_input,const VECTOR<T,3>& color2_input,const int octaves,const T lacunarity,const T gain,const T low,const T high,const T vain_value,
        const T vain_width,const T clamp_width,RENDER_WORLD<T>& world) 
        :MATERIAL_SHADER<T>(world),color1(color1_input),color2(color2_input),octaves(octaves),lacunarity(lacunarity),gain(gain),low(low),high(high),vain_value(vain_value),vain_width(vain_width),
        clamp_width(clamp_width)
    {
        RANGE<VECTOR<T,1> > dark(vain_value-vain_width,vain_value+vain_width);
        T clamp_value=vain_value+clamp_width;
        curve.Add_Control_Point(0,clamp_value);
        curve.Add_Control_Point(dark.min_corner.x-clamp_width,clamp_value);
        curve.Add_Control_Point(dark.min_corner.x,dark.max_corner.x);
        curve.Add_Control_Point(dark.Size()/4+dark.min_corner.x,(T).3);
        curve.Add_Control_Point(dark.Size()/2+dark.min_corner.x,(T).2);
        curve.Add_Control_Point(3*dark.Size()/4+dark.min_corner.x,(T).3);
        curve.Add_Control_Point(dark.max_corner.x,dark.max_corner.x);
        curve.Add_Control_Point(dark.max_corner.x+clamp_width,clamp_value);
    }

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {RAY<VECTOR<T,3> > object_space_ray=intersection_object.Object_Space_Ray(ray.ray);VECTOR<T,3> object_space_position=object_space_ray.Point(object_space_ray.t_max);
    TV p=intersection_object.Get_Solid_Texture_Coordinates(object_space_position,ray.ray.aggregate_id);
    T a=curve.Value(Smooth((T)NOISE<T>::Turbulence(p,octaves,lacunarity,gain),low,high)); 
    return (1-a)*color1+a*color2;}

//#####################################################################
};
}
#endif

