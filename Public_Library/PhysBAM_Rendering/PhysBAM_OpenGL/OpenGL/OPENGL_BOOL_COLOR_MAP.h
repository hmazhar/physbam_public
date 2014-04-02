//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BOOL_COLOR_MAP
//#####################################################################
#ifndef __OPENGL_BOOL_COLOR_MAP__
#define __OPENGL_BOOL_COLOR_MAP__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
namespace PhysBAM{

class OPENGL_BOOL_COLOR_MAP : public OPENGL_COLOR_MAP<bool>
{
public:
    OPENGL_COLOR true_color;
    OPENGL_COLOR false_color;

    OPENGL_BOOL_COLOR_MAP(const OPENGL_COLOR &true_color,const OPENGL_COLOR &false_color)
        :true_color(true_color),false_color(false_color)
    {}

    virtual OPENGL_COLOR Lookup(bool x) const PHYSBAM_OVERRIDE
    {return x?true_color:false_color;}
};
}
#endif
