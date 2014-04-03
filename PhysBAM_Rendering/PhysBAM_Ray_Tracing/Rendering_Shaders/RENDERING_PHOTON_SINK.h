//#####################################################################
// Copyright 2003-2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_PHOTON_SINK__
#define __RENDERING_PHOTON_SINK__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_PHOTON_SINK:public MATERIAL_SHADER<T>
{
public:
    using MATERIAL_SHADER<T>::world;

    RENDERING_PHOTON_SINK(RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input)
    {}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {if(VECTOR<T,3>::Dot_Product(same_side_normal,intersection_object.Normal(intersection_point,ray.ray.aggregate_id))>0){
        RENDERING_RAY<T> render_ray=RENDERING_RAY<T>(RAY<VECTOR<T,3> >(intersection_point-same_side_normal*intersection_object.small_number*4,ray.ray.direction),ray.ray_contribution,&entering_object);
        return world.Cast_Ray(render_ray,ray);}
    else return world.background_color;}

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
        unsigned int photon_flags) const
    {return;}

//#####################################################################
};
}
#endif
