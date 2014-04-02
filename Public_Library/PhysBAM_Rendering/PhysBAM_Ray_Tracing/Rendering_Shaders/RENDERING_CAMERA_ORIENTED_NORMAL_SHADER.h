//#####################################################################
// Copyright 2007, Justin Solomon.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_CAMERA_ORIENTED_NORMAL_SHADER__
#define __RENDERING_CAMERA_ORIENTED_NORMAL_SHADER__
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_SEGMENTED_CURVE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_CAMERA_ORIENTED_NORMAL_SHADER:public MATERIAL_SHADER<T>
{
public:
    using MATERIAL_SHADER<T>::world;
    const MATERIAL_SHADER<T>& shader;

    RENDERING_CAMERA_ORIENTED_NORMAL_SHADER(const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),shader(shader_input)
    {}

    VECTOR<T,3> Rotate_Normal_Toward_Camera(const int aggregate_id,const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3> &intersection_point, const VECTOR<T,3> &normal) const
    {if(const RENDERING_SEGMENTED_CURVE<T> *rendering_segmented_curve=dynamic_cast<const RENDERING_SEGMENTED_CURVE<T>*>(&intersection_object)){
        const VECTOR<int,2>& nodes=rendering_segmented_curve->segmented_curve.mesh.elements(aggregate_id);
        VECTOR<T,3> hair_direction=rendering_segmented_curve->World_Space_Vector(rendering_segmented_curve->segmented_curve.particles.X(nodes[2])-rendering_segmented_curve->segmented_curve.particles.X(nodes[1]));
        return (world.camera.position-intersection_point).Projected_Orthogonal_To_Unit_Direction(hair_direction.Normalized()).Normalized();}
    else{LOG::cout<<"camera-oriented normal shader only can be used on segmented curves"<<std::endl;PHYSBAM_FATAL_ERROR();}}
    
    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {return shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,
                Rotate_Normal_Toward_Camera(ray.ray.aggregate_id,intersection_object,intersection_point,same_side_normal));}

    VECTOR<T,3> Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {return shader.Shade_Surface_Using_Indirect_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,
                Rotate_Normal_Toward_Camera(ray.ray.aggregate_id,intersection_object,intersection_point,same_side_normal));}

    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,
        const RENDERING_RAY<T>& full_ray) const
    {return shader.Shade_Light_Ray(ray,exiting_object,entering_object,intersection_object,intersection_point,
                Rotate_Normal_Toward_Camera(ray.ray.aggregate_id,intersection_object,intersection_point,same_side_normal),light,full_ray);}

    VECTOR<T,3> Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {return shader.Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,
                Rotate_Normal_Toward_Camera(ray.ray.aggregate_id,intersection_object,intersection_point,same_side_normal));}

    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const 
    {return shader.Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,
            Rotate_Normal_Toward_Camera(ray.ray.aggregate_id,intersection_object,intersection_point,same_side_normal),power,type,diffuse_bounces,specular_bounces);}
//#####################################################################
};
}
#endif
