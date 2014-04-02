//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW_PBUFFER.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#ifdef WIN32
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WGL_PBUFFER.h>
#elif defined(__APPLE__)
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AGL_PBUFFER.h>
#else
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GLX_PBUFFER.h>
#endif

using namespace PhysBAM;

//#####################################################################
// OPENGL_WINDOW_PBUFFER
//#####################################################################
OPENGL_WINDOW_PBUFFER::
OPENGL_WINDOW_PBUFFER(OPENGL_WORLD& opengl_world_input,const std::string& window_title_input,const int width_input,const int height_input)
    :OPENGL_WINDOW(opengl_world_input),width(width_input),height(height_input)
{
    pbuffer=new OPENGL_PBUFFER;
    if(!pbuffer->Create(width,height)) PHYSBAM_FATAL_ERROR("Could not set up pbuffer");
    // need glut init to make sure we can run glut functions
    static int argc=1;static const char *(argv[1]);argv[0]="Visualization";
    glutInit(&argc,(char**)argv);
}
//#####################################################################
// ~OPENGL_WINDOW_PBUFFER
//#####################################################################
OPENGL_WINDOW_PBUFFER::
~OPENGL_WINDOW_PBUFFER()
{
    delete pbuffer;
}
//#####################################################################
// Function Handle_Idle
//#####################################################################
void OPENGL_WINDOW_PBUFFER::
Setup_Idle(const bool use)
{}
//#####################################################################
// Function Handle_Idle
//#####################################################################
void OPENGL_WINDOW_PBUFFER::
Setup_Timer(const float wait_milliseconds)
{}
//#####################################################################
// Function Handle_Idle
//#####################################################################
void OPENGL_WINDOW_PBUFFER::
Redisplay()
{}
//#####################################################################
// Function Handle_Idle
//#####################################################################
void OPENGL_WINDOW_PBUFFER::
Main_Loop()
{}
//#####################################################################
// Request_Resize
//#####################################################################
void OPENGL_WINDOW_PBUFFER::
Request_Resize(const int width,const int height)
{}
//#####################################################################
// Request_Move
//#####################################################################
void OPENGL_WINDOW_PBUFFER::
Request_Move(const int x,const int y)
{}
//#####################################################################
// Width
//#####################################################################
int OPENGL_WINDOW_PBUFFER::
Width() const
{ 
    return width;
}
//#####################################################################
// Height
//#####################################################################
int OPENGL_WINDOW_PBUFFER::
Height() const
{
    return height;
}
//#####################################################################
