//#####################################################################
// Copyright 2003, 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_LIGHT_SHADER
//#####################################################################
#ifndef __RENDERING_LIGHT_SHADER__
#define __RENDERING_LIGHT_SHADER__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_LIGHT_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
protected:
    RENDERING_LIGHT<T>& light;

public:
    RENDERING_LIGHT_SHADER(RENDERING_LIGHT<T>& light_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),light(light_input)
    {}

    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
    {return light.Emitted_Light(ray)*(4*(T)pi);} // we want outgoing radiance here, not the irradiance at a surface

    TV Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal,const RENDERING_LIGHT<T>& light,
        const RENDERING_RAY<T>& full_ray) const 
    {return light.Emitted_Light(ray);} // we want the irradiance due to the light source in this case

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal,const TV& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
    {return;}
    
    virtual TV Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const 
    {return TV();}

//#####################################################################
};
}
#endif
