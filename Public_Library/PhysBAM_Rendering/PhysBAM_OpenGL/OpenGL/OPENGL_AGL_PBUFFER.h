//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_AGL_PBUFFER
//#####################################################################
#ifndef __OPENGL_AGL_PBUFFER__
#define __OPENGL_AGL_PBUFFER__

#ifdef __APPLE__
#include <AGL/agl.h>
#define OPENGL_AGL_PBUFFER_SUPPORTED
#endif

#ifdef OPENGL_AGL_PBUFFER_SUPPORTED
namespace PhysBAM{
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
        AGLPbuffer pbuffer;        
    };
}
#else
#include <PhysBAM_Tools/Log/LOG.h>
namespace PhysBAM {
    class OPENGL_PBUFFER
    {
    public:
        bool Create(int width, int height) { LOG::cerr << "OPENGL_AGL_PBUFFER not supported" << std::endl; return false; }
        void Destroy() { }
    };
}
#endif
#endif
