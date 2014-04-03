//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// pbuffers code taken from pbdemo.c and pbutil.{h,c} in
// http://cvs.sourceforge.net/viewcvs.py/mesa3d/Mesa/xdemos/
// written by Brain Paul (April 1997; Updated on 5 October 2002)
//#####################################################################
#ifdef __APPLE__
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AGL_PBUFFER.h>
#include <cstdio>
#include <cstring>
#include <iostream>

#ifdef OPENGL_AGL_PBUFFER_SUPPORTED

namespace PhysBAM {

OPENGL_PBUFFER::OPENGL_PBUFFER()
:verbose(false),pbuffer(0)
{}

OPENGL_PBUFFER::~OPENGL_PBUFFER()
{
    Destroy();
}

bool OPENGL_PBUFFER::
Create(int width, int height)
{
    if (pbuffer) {
        LOG::cerr << "Destroying old pbuffer before creating new one" << std::endl;
        Destroy();
    }

    GLint attributes[] = {AGL_RGBA,AGL_DOUBLEBUFFER,AGL_NONE};
    AGLPixelFormat format = aglChoosePixelFormat(NULL, 0, attributes);
    if(format==NULL){
        LOG::cout<<"Could not set up the pixel format"<<std::endl;return false;}
    
    AGLContext context = aglCreateContext(format, NULL);
    if(context==NULL){
        LOG::cout<<"Could not set up the pbuffer context"<<std::endl;return false;}
    aglSetCurrentContext(context);

    if(!aglCreatePBuffer(width, height, GL_TEXTURE_2D, GL_RGBA, 0, &pbuffer)){
        LOG::cout<<"Error: couldn't create pbuffer"<<std::endl;
        LOG::cout<<aglErrorString(aglGetError());
        return false;
    }

    if(!aglSetPBuffer(context, pbuffer, 0, 0, aglGetVirtualScreen(context))){
        LOG::cout<<"Error: couldn't set pbuffer"<<std::endl;
        LOG::cout<<aglErrorString(aglGetError());
        return false;
    }
       
    return true;  /* Success!! */
}

void OPENGL_PBUFFER::
Destroy()
{
    if(pbuffer){
        if(verbose) LOG::cout<<"Destroying pbuffer"<<std::endl;
        aglDestroyPBuffer(pbuffer);}
}
}

#endif
#endif
