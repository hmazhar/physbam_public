//#####################################################################
// Copyright 2004, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_BUMP_MAP_SHADER__
#define __RENDERING_BUMP_MAP_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BUMP_MAP_SHADER:public MATERIAL_SHADER<T>
{
public:
    const MATERIAL_SHADER<T>* shader;
    
    RENDERING_BUMP_MAP_SHADER(const MATERIAL_SHADER<T>* shader_input,RENDER_WORLD<T>& world_input) 
        :MATERIAL_SHADER<T>(world_input),shader(shader_input)
    {}
    
    virtual VECTOR<T,3> Perturbed_Normal( const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& object_space_position,const T s,const T t,const VECTOR<T,3>& object_space_normal,const VECTOR<T,3>& object_space_tangent,const VECTOR<T,3>& object_space_bitangent) const=0;
    
    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {VECTOR<T,3> object_space_intersection=intersection_object.Object_Space_Point(intersection_point);
    VECTOR<T,3> object_space_normal=intersection_object.Object_Space_Vector(same_side_normal);
    T s,t;intersection_object.Get_Texture_Coordinates(object_space_intersection,ray.ray.aggregate_id,s,t);
    VECTOR<T,3> tangent,bitangent;intersection_object.Get_Object_Space_Tangent_And_Bitangent(object_space_intersection,object_space_normal,ray.ray.aggregate_id,tangent,bitangent);
    VECTOR<T,3> perturbed_normal=intersection_object.World_Space_Vector(Perturbed_Normal(intersection_object,object_space_intersection,s,t,object_space_normal,tangent,bitangent));
    return shader->Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,perturbed_normal);}

//#####################################################################
};
}
#endif
