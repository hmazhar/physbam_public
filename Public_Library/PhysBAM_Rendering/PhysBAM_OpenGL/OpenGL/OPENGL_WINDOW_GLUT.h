//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_WINDOW_GLUT__
#define __OPENGL_WINDOW_GLUT__

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW.h>
namespace PhysBAM{

class OPENGL_WINDOW_GLUT:public OPENGL_WINDOW
{
    using OPENGL_WINDOW::opengl_world;
    static OPENGL_WINDOW_GLUT* single;
    int main_window;
    int width,height;

//#####################################################################
public:
    OPENGL_WINDOW_GLUT(OPENGL_WORLD& world_input,const std::string& window_title_input,const int width_input,const int height_input);
    virtual ~OPENGL_WINDOW_GLUT();
    void Setup_Idle(const bool use) PHYSBAM_OVERRIDE;
    void Setup_Timer(const float wait_milliseconds) PHYSBAM_OVERRIDE;
    void Redisplay() PHYSBAM_OVERRIDE;
    void Main_Loop() PHYSBAM_OVERRIDE;
    void Request_Resize(const int width,const int height) PHYSBAM_OVERRIDE;
    void Request_Move(const int x,const int y) PHYSBAM_OVERRIDE;
    int Width() const PHYSBAM_OVERRIDE;
    int Height() const PHYSBAM_OVERRIDE;
private:
    static void Handle_Idle_Glut();
    static void Handle_Timer_Glut(int value);
    static void Handle_Display_Glut();
    static void Handle_Reshape_Glut(int w,int h);
    static void Handle_Special_Keypress_Glut(int key,int x,int y);
    static void Handle_Keypress_Glut(unsigned char key,int x,int y);
    static void Handle_Click_Glut(int button,int state,int x,int y);
    static void Handle_Drag_Glut(int x,int y);
//#####################################################################
};
}
#endif
