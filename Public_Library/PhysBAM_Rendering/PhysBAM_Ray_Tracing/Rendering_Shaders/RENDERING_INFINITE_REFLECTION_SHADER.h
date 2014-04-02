//#####################################################################
// Copyright 2004-2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_INFINITE_REFLECTION_SHADER  
//#####################################################################
#ifndef __RENDERING_INFINITE_REFLECTION_SHADER__
#define __RENDERING_INFINITE_REFLECTION_SHADER__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Math_Tools/wrap.h>
#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_INFINITE_REFLECTION_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<T,2> TV2;
public:
    GRID<TV2> grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > pixels;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> > interpolation;
    T rotation; // stored as a value between -pi and pi
    
    RENDERING_INFINITE_REFLECTION_SHADER(RENDER_WORLD<T>& world_input,const T rotation_input) 
        :MATERIAL_SHADER<T>(world_input)
    {rotation=wrap(rotation_input,(T)-pi,(T)pi);}

    void Initialize(const std::string& filename,const T scale=1)
    {IMAGE<T>::Read(filename,pixels);grid.Initialize(pixels.counts,RANGE<VECTOR<T,2> >::Unit_Box());pixels*=scale;}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {static T one_over_two_pi=(T)1/((T)two_pi),one_over_pi=(T)1/(T)pi;
    T theta=acos(ray.ray.direction.y),t=1-theta*one_over_pi; // maps theta=0 to t=1, theta=pi to t=0
    T phi=atan2(ray.ray.direction.x,ray.ray.direction.z)-rotation,s=phi*one_over_two_pi;if(s<0)s+=1; // maps phi=0 to s=0, s increases with positive rotation about y axis
    return interpolation.Clamped_To_Array_Cell(grid,pixels,VECTOR<T,2>(s,t));}
    
//#####################################################################
};
}
#endif

