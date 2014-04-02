//#####################################################################
// Copyright 2004, Frank Losasso, Andrew Selle, Jiayi Chong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_REFLECTION_MAP_SHADER__
#define __RENDERING_REFLECTION_MAP_SHADER__

#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_TEXTURE_IMAGE_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_REFLECTION_MAP_SHADER:public RENDERING_TEXTURE_IMAGE_SHADER<T>
{
public:
    using RENDERING_TEXTURE_IMAGE_SHADER<T>::world;

    RENDERING_REFLECTION_MAP_SHADER(RENDER_WORLD<T>& world_input) 
        :RENDERING_TEXTURE_IMAGE_SHADER<T>(world_input)
    {}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
    RENDERING_RAY<T> outgoing_ray(RAY<VECTOR<T,3> >(same_side_position,world.Reflected_Direction(same_side_normal,ray.ray.direction),true),1,&entering_object);
    return world.Cast_Ray(outgoing_ray,ray)*RENDERING_TEXTURE_IMAGE_SHADER<T>::Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,
        intersection_point,same_side_normal);}

//#####################################################################
};
}
#endif
