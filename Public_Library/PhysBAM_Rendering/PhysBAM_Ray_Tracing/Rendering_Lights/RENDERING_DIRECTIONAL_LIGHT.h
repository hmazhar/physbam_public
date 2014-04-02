//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_DIRECTIONAL_LIGHT__
#define __RENDERING_DIRECTIONAL_LIGHT__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
namespace PhysBAM{

template<class T>
class RENDERING_DIRECTIONAL_LIGHT:public RENDERING_LIGHT<T>
{
public:
    using RENDERING_LIGHT<T>::world;using RENDERING_LIGHT<T>::position;using RENDERING_LIGHT<T>::color;using RENDERING_LIGHT<T>::brightness;

    VECTOR<T,3> direction;

    RENDERING_DIRECTIONAL_LIGHT(const VECTOR<T,3>& direction_input,const VECTOR<T,3>& color_input,const T brightness_input,RENDER_WORLD<T>& world_input) 
        :RENDERING_LIGHT<T>(direction_input,color_input,brightness_input,world_input),direction(direction_input)
    {
        direction.Normalize();
    }

    void Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,ARRAY<RAY<VECTOR<T,3> > >& sample_array)const PHYSBAM_OVERRIDE
    {sample_array.Resize(1);sample_array(1)=RAY<VECTOR<T,3> >(surface_position,-direction,true);}

    VECTOR<T,3> Emitted_Light(const RENDERING_RAY<T>& ray)const PHYSBAM_OVERRIDE
    {return brightness*color;}

//#####################################################################
};
}
#endif
