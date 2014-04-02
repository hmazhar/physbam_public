//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// class OPENGL_COLOR_RAMP
//##################################################################### 
#ifndef __OPENGL_COLOR_RAMP__
#define __OPENGL_COLOR_RAMP__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
namespace PhysBAM{

template<class T>
class OPENGL_COLOR_RAMP : public OPENGL_COLOR_MAP<T>
{
public:

    ARRAY<T> color_x;
    ARRAY<OPENGL_COLOR> less_colors,equal_colors,greater_colors;

    OPENGL_COLOR_RAMP();
    ~OPENGL_COLOR_RAMP();

    void Add_Color(T x,const OPENGL_COLOR& color)
    {Add_Color(x,color,color,color);}

    void Add_Color(T x,const OPENGL_COLOR& color,const OPENGL_COLOR& exact_color)
    {Add_Color(x,color,exact_color,color);}

    OPENGL_COLOR Lookup(T x) const PHYSBAM_OVERRIDE;
    void Add_Color(T x,const OPENGL_COLOR& less_color,const OPENGL_COLOR& exact_color,const OPENGL_COLOR& greater_color);
    static OPENGL_COLOR_RAMP<T> *Matlab_Jet(T value_min,T value_max);
    static OPENGL_COLOR_RAMP<T> *Matlab_Hot(T value_min,T value_max);
    static OPENGL_COLOR_RAMP<T> *Two_Color_Ramp(T value_min,T value_max,const OPENGL_COLOR& color_min,const OPENGL_COLOR& color_min_exact,
        const OPENGL_COLOR& color_max,const OPENGL_COLOR& color_max_exact);
    static OPENGL_COLOR_RAMP<T> *Two_Color_Ramp(T value_min,T value_max,const OPENGL_COLOR& color_min,const OPENGL_COLOR& color_max);
    static OPENGL_COLOR_RAMP<T> *Levelset_Color_Constant_Ramp(const OPENGL_COLOR& negative_color,const OPENGL_COLOR& positive_color);
    static OPENGL_COLOR_RAMP<T> *Levelset_Color_Linear_Ramp(const OPENGL_COLOR& negative_color,const OPENGL_COLOR& positive_color,T abs_value_max);
//#####################################################################
};   
}
#endif
