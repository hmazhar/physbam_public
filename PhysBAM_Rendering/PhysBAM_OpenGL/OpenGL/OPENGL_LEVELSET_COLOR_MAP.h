//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_LEVELSET_COLOR_MAP
//#####################################################################
#ifndef __OPENGL_LEVELSET_COLOR_MAP__
#define __OPENGL_LEVELSET_COLOR_MAP__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
namespace PhysBAM{

template<class T>
class OPENGL_LEVELSET_COLOR_MAP : public OPENGL_COLOR_MAP<T>
{
public:
    OPENGL_COLOR negative_color;
    OPENGL_COLOR positive_color;

    OPENGL_LEVELSET_COLOR_MAP(const OPENGL_COLOR& negative_color_input,const OPENGL_COLOR& positive_color_input)
        :negative_color(negative_color_input),positive_color(positive_color_input)
    {}

    virtual OPENGL_COLOR Lookup(T x) const PHYSBAM_OVERRIDE
    {if(x<=0) return (1+2*x)*negative_color;else return (1-2*x)*positive_color;}
};
}
#endif
