//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_TEXTURE_IMAGE_SHADER  
//#####################################################################
#ifndef __RENDERING_TEXTURE_IMAGE_SHADER__
#define __RENDERING_TEXTURE_IMAGE_SHADER__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_TEXTURE_IMAGE_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    GRID<VECTOR<T,2> > grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > pixels;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> > linear_interpolation;
    CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> > cubic_interpolation;
    INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> >* interpolation;
    bool wrap_s,wrap_t;
    
    RENDERING_TEXTURE_IMAGE_SHADER(RENDER_WORLD<T>& world_input,const bool wrap_s_input=true,const bool wrap_t_input=true,const bool cubic_interpolation_input=false)
        :MATERIAL_SHADER<T>(world_input),wrap_s(wrap_s_input),wrap_t(wrap_t_input)
    {
        if(cubic_interpolation_input) interpolation=&cubic_interpolation;else interpolation=&linear_interpolation;
    }

    bool Valid() const
    {return grid.counts.x>0 && grid.counts.y>0 && grid.counts==pixels.counts;}

//#####################################################################
    void Setup_Interpolation(const int m,const int n);
    void Initialize(const std::string& filename);
    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
         const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const;
//#####################################################################
};
}
#endif

