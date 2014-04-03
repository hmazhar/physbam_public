//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_WINDOW__
#define __OPENGL_WINDOW__
#include <string>
namespace PhysBAM{

class OPENGL_WORLD;

class OPENGL_WINDOW
{
public:
    OPENGL_WORLD& opengl_world;
    std::string window_title;

    OPENGL_WINDOW(OPENGL_WORLD& opengl_world_input)
        :opengl_world(opengl_world_input),window_title("PhysBAM OpenGL")
    {}

    virtual ~OPENGL_WINDOW()
    {}

//#####################################################################
    virtual void Setup_Idle(const bool use)=0;
    virtual void Setup_Timer(const float wait_milliseconds)=0;
    virtual void Redisplay()=0;
    virtual void Main_Loop()=0;
    virtual void Request_Resize(const int width,const int height)=0;
    virtual void Request_Move(const int x,const int y)=0;
    virtual int Width() const=0;
    virtual int Height() const=0;
//#####################################################################
};
}
#endif
