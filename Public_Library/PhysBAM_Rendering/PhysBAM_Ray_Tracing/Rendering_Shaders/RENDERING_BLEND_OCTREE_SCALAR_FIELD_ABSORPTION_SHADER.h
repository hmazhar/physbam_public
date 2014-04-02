#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2003, 2004, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_BLEND_OCTREE_SCALAR_FIELD_ABSORPTION_SHADER__
#define __RENDERING_BLEND_OCTREE_SCALAR_FIELD_ABSORPTION_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OCTREE_IMPLICIT_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BLEND_OCTREE_SCALAR_FIELD_ABSORPTION_SHADER:public VOLUMETRIC_SHADER<T>
{
public:
    T absorption_coefficient1,absorption_coefficient2;
    const VECTOR<T,3> absorption_spectrum1,absorption_spectrum2;

    const ARRAY<T>* values1;
    const ARRAY<T>* values2; // may be null
    bool use_abs_value;
    T tolerance;
    T blend_band;
    int blend_mode;

    RENDERING_BLEND_OCTREE_SCALAR_FIELD_ABSORPTION_SHADER(const T absorption_coefficient1_input,const T absorption_coefficient2_input,
                                                          const VECTOR<T,3>& absorption_spectrum1_input,const VECTOR<T,3>& absorption_spectrum2_input,
                                                          const ARRAY<T>* values1_input,const ARRAY<T>* values2_input,RENDER_WORLD<T>& world_input) 
        :VOLUMETRIC_SHADER<T>(world_input),absorption_coefficient1(absorption_coefficient1_input),absorption_coefficient2(absorption_coefficient2_input),absorption_spectrum1(absorption_spectrum1_input),absorption_spectrum2(absorption_spectrum2_input),
        values1(values1_input),values2(values2_input),use_abs_value(false),tolerance((T)1e-8),blend_band((T).1),blend_mode(1)
    {}

    VECTOR<T,3> Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const VECTOR<T,3>& color) PHYSBAM_OVERRIDE
    {VECTOR<T,3> sample_point=ray.ray.Point(ray.ray.t_max/2);
    if(ray.ray.semi_infinite||!object.Inside(sample_point)) return color;
    const OCTREE_GRID<T>& grid=((const RENDERING_OCTREE_IMPLICIT_SURFACE<T>*)&object)->implicit_surface->levelset.grid;
    const OCTREE_CELL<T>* cell=grid.Leaf_Cell(ray.ray.Point(ray.ray.t_max/2));
    T value1=LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,cell,*values1,sample_point);
    T blend=0;
    if(blend_mode==1){assert(values2);
        T value2=LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,cell,*values2,sample_point);
        if(use_abs_value){value1=abs(value1);value2=abs(value2);}
        blend=clamp((T)0.5+(value1-value2)/blend_band,(T)0,(T)1);}
    else{assert(blend_mode==2);
        blend=clamp(value1/blend_band,(T)0,(T)1);}

    T absorption_coefficient=(1-blend)*absorption_coefficient1+blend*absorption_coefficient2;
    VECTOR<T,3> absorption_spectrum=(1-blend)*absorption_spectrum1+blend*absorption_spectrum2;

    VECTOR<T,3> attenuation=exp(-ray.ray.t_max*absorption_coefficient*absorption_spectrum);
    if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("attenuation=%f %f %f\n",attenuation.x,attenuation.y,attenuation.z));
    return color*attenuation;}

    VECTOR<T,3> Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,const VECTOR<T,3>& light_color) PHYSBAM_OVERRIDE
    {return Attenuate_Color(ray,object,light_color);}

//#####################################################################
};
}
#endif
#endif
