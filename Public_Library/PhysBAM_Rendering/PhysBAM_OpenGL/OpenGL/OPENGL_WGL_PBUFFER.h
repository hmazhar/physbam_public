//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_WGL_PBUFFER
//#####################################################################
#ifndef __OPENGL_WGL_PBUFFER__
#define __OPENGL_WGL_PBUFFER__

#if defined(USERNAME_aselle)
#define OPENGL_WGL_PBUFFER_SUPPORTED
#endif
#if defined(WIN32) && defined(OPENGL_WGL_PBUFFER_SUPPORTED) 

#include <PhysBAM_Tools/Log/LOG.h>
#include <GL/gl.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/glATI.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/wglATI.h>
#include <windows.h>

namespace PhysBAM
{
    class OPENGL_PBUFFER
    {
    public:
        OPENGL_PBUFFER()
        {}

        ~OPENGL_PBUFFER()
        {}

        bool Create(int width, int height);
        void Destroy();

    public:
        bool verbose;

    private:
        // stupid window
        HWND hwnd;
        // device and render context for frame buffer
        HDC hGLDC;
        HGLRC hGLRC,hPBufferGLRC;

        // handle, device context and render context for pbuffer buffer
        HPBUFFERARB hPBuffer;
        HDC hPBufferDC;

        // variables to hold pbuffers dimensions
        int pbufferWidth;
        int pbufferHeight;

        // texture object ID for binding pbuffer to a texture
        GLuint pbufferTextureID;

        // allows us to foolproof the interface to prevent using the pbuffer before
        // it is initialized
        static bool pbufferInitialized;


    };
}
#else
#include <PhysBAM_Tools/Log/LOG.h>
#include <iostream>
namespace PhysBAM
{
    class OPENGL_PBUFFER
    {
    public:
        bool Create(int width, int height) { LOG::cerr << "OPENGL_WGL_PBUFFER not supported" << std::endl; return false; }
        void Destroy() { }
    };
}


#endif

#endif // _opengl_wgl_pbuffer_h_
