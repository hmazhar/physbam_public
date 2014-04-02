//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GLX_PBUFFER
//#####################################################################
#ifndef __OPENGL_GLX_PBUFFER__
#define __OPENGL_GLX_PBUFFER__

#ifdef __linux__
#ifndef GLX_GLXEXT_PROTOTYPES
#define GLX_GLXEXT_PROTOTYPES
#endif
#include <GL/glx.h>
//#if defined(GLX_VERSION_1_3) && defined(GLX_SGIX_pbuffer)
#if defined(GLX_VERSION_1_3)
#define OPENGL_GLX_PBUFFER_SUPPORTED
#endif
#endif
#ifndef OPENGL_GLX_PBUFFER_SUPPORTED
#include <PhysBAM_Tools/Log/LOG.h>
#endif

namespace PhysBAM{

#ifdef OPENGL_GLX_PBUFFER_SUPPORTED
class OPENGL_PBUFFER
{
public:
    OPENGL_PBUFFER();
    ~OPENGL_PBUFFER();

    bool Create(int width, int height);
    void Destroy();

public:
    bool verbose;

private:
    Display *display;
    GLXFBConfig fbconfig;
    GLXPbuffer pbuffer;


};
#else

class OPENGL_PBUFFER
{
public:
    bool Create(int width, int height) { LOG::cerr << "OPENGL_GLX_PBUFFER not supported" << std::endl; return false; }
    void Destroy() { }
};
#endif

}

#endif
