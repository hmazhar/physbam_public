//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TEXTURE
//#####################################################################
#ifndef __OPENGL_TEXTURE__
#define __OPENGL_TEXTURE__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>

namespace PhysBAM
{
class OPENGL_COLOR;

class OPENGL_TEXTURE
{
public:
    GLuint  id;
    int     padded_width, padded_height;    // to make it a power of 2
    int     width, height;
    double  min_s, min_t, max_s, max_t;
    bool    smooth_shading;

    OPENGL_TEXTURE();
    ~OPENGL_TEXTURE();

    void Initialize(int width_input, int height_input, double border=0);
    void Update_Texture(const OPENGL_COLOR *bitmap);
    void Set_Smooth_Shading(bool smooth_shading_input=true);
    void Toggle_Smooth_Shading();
};

}

#endif
