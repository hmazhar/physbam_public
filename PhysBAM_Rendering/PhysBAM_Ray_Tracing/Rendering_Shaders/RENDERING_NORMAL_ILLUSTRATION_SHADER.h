//#####################################################################
// Copyright 2003, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_NORMAL_ILLUSTRATION_SHADER__
#define __RENDERING_NORMAL_ILLUSTRATION_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_NORMAL_ILLUSTRATION_SHADER:public MATERIAL_SHADER<T>
{
public:
    RENDERING_NORMAL_ILLUSTRATION_SHADER(RENDER_WORLD<T>& world) 
        :MATERIAL_SHADER<T>(world)
    {}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {return (same_side_normal+VECTOR<T,3>(1,1,1))*(T).5;}

//#####################################################################
};
}
#endif
