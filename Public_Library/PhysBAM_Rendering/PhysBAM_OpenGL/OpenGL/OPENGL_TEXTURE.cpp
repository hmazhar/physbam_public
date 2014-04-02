//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TEXTURE.h>
#include <cmath>
#include <cstring>
#include <iostream>
using namespace PhysBAM;
using namespace std;

namespace
{
    int Round_To_Power_Of_Two(int n)
    {
        return 1 << (int)ceil(log((double)n)/log((double)2.0));
    }
}

OPENGL_TEXTURE::OPENGL_TEXTURE()
{
    glGenTextures(1, &id);
}

OPENGL_TEXTURE::~OPENGL_TEXTURE()
{
    glDeleteTextures(1, &id);
}

void
OPENGL_TEXTURE::Initialize(int width_input, int height_input, double border)
{
    width = width_input;
    height = height_input;

    padded_width = Round_To_Power_Of_Two(width);
    padded_height = Round_To_Power_Of_Two(height);

    min_s = border/padded_width;
    max_s = (double)(width-border)/padded_width;
    min_t = border/padded_width;
    max_t = (double)(height-border)/padded_height;

#ifdef DEBUG
    cout << "Texture: got actual " << width << "x" << height 
         << ", padded " << padded_width << "x" << padded_height 
         << " maxs " << max_s << " maxt " << max_t << endl;
#endif

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glBindTexture(GL_TEXTURE_2D, id);

    Set_Smooth_Shading(false);

    OPENGL_COLOR *initial_bitmap = new OPENGL_COLOR[padded_width*padded_height];
    memset(initial_bitmap, 0, padded_width*padded_height*sizeof(OPENGL_COLOR));
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, padded_width, padded_height, 0, GL_RGBA, GL_FLOAT, initial_bitmap);
    delete[] initial_bitmap;
}

void
OPENGL_TEXTURE::Update_Texture(const OPENGL_COLOR *bitmap)
{
    glBindTexture(GL_TEXTURE_2D, id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_FLOAT, bitmap);
}

void
OPENGL_TEXTURE::Set_Smooth_Shading(bool smooth_shading_input)
{
    smooth_shading = smooth_shading_input;
    glBindTexture(GL_TEXTURE_2D,id);
    if(smooth_shading) {    
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }
}

void
OPENGL_TEXTURE::Toggle_Smooth_Shading()
{
    Set_Smooth_Shading(!smooth_shading);
}
