//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_CONSTANT_COLOR_MAP
//#####################################################################
#ifndef __OPENGL_CONSTANT_COLOR_MAP__
#define __OPENGL_CONSTANT_COLOR_MAP__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
namespace PhysBAM{

template<class T>
class OPENGL_CONSTANT_COLOR_MAP : public OPENGL_COLOR_MAP<T>
{
public:
    OPENGL_COLOR color;

    OPENGL_CONSTANT_COLOR_MAP(const OPENGL_COLOR &color)
        :color(color)
    {}

    virtual OPENGL_COLOR Lookup(T x) const PHYSBAM_OVERRIDE
    {return color;}
};
}
#endif
