//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COLOR_MAP
//#####################################################################
#ifndef __OPENGL_COLOR_MAP__
#define __OPENGL_COLOR_MAP__

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
namespace PhysBAM{

template<class T>
class OPENGL_COLOR_MAP
{
public:
    virtual ~OPENGL_COLOR_MAP(){}

    virtual OPENGL_COLOR Lookup(T x) const=0;
};
}
#endif
