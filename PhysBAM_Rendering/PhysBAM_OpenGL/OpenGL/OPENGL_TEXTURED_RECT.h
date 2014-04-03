//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TEXTURED_RECT
//#####################################################################
#ifndef __OPENGL_TEXTURED_RECT__
#define __OPENGL_TEXTURED_RECT__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TEXTURE.h>

namespace PhysBAM
{

class OPENGL_TEXTURED_RECT : public OPENGL_OBJECT
{
public:
    double width, height;
    OPENGL_TEXTURE *texture;

    OPENGL_TEXTURED_RECT();

    void Set_Texture(OPENGL_TEXTURE *texture_input);
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
