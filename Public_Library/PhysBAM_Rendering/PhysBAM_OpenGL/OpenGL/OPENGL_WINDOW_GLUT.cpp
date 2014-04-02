//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW_GLUT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>

using namespace PhysBAM;

OPENGL_WINDOW_GLUT* OPENGL_WINDOW_GLUT::single=0;

//#####################################################################
// OPENGL_WINDOW_GLUT
//#####################################################################
OPENGL_WINDOW_GLUT::
OPENGL_WINDOW_GLUT(OPENGL_WORLD& opengl_world_input,const std::string& window_title_input,const int width_input,const int height_input)
    :OPENGL_WINDOW(opengl_world_input),width(width_input),height(height_input)
{
    if(single) PHYSBAM_FATAL_ERROR("Only one glut context allowed");
    single=this;

    static int argc=1;static const char *(argv[1]);argv[0]="Visualization";
    glutInit(&argc,(char**)argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_ALPHA);
    glutInitWindowSize(width,height);
    main_window=glutCreateWindow(window_title_input.c_str());
    glutDisplayFunc(Handle_Display_Glut);
    glutReshapeFunc(Handle_Reshape_Glut);
    glutKeyboardFunc(Handle_Keypress_Glut);
    glutSpecialFunc(Handle_Special_Keypress_Glut);
    glutMouseFunc(Handle_Click_Glut);
    glutMotionFunc(Handle_Drag_Glut);
}
//#####################################################################
// ~OPENGL_WINDOW_GLUT
//#####################################################################
OPENGL_WINDOW_GLUT::
~OPENGL_WINDOW_GLUT()
{single=0;}

//#####################################################################
// Setup_Idle
//#####################################################################
void OPENGL_WINDOW_GLUT::
Setup_Idle(const bool use)
{
    glutIdleFunc(use?Handle_Idle_Glut:0);
}
//#####################################################################
// Setup_Timer
//#####################################################################
void OPENGL_WINDOW_GLUT::
Setup_Timer(const float wait_milliseconds)
{
    glutTimerFunc((int)(wait_milliseconds*1000)+1,Handle_Timer_Glut,0);
}
//#####################################################################
// Redisplay
//#####################################################################
void OPENGL_WINDOW_GLUT::
Redisplay()
{
    glutPostRedisplay();
}
//#####################################################################
// Main_Loop
//#####################################################################
void OPENGL_WINDOW_GLUT::
Main_Loop()
{
    glutMainLoop();
}
//#####################################################################
// Request_Resize
//#####################################################################
void OPENGL_WINDOW_GLUT::
Request_Resize(const int width,const int height)
{
    glutReshapeWindow(width,height);
}
//#####################################################################
// Request_Move
//#####################################################################
void OPENGL_WINDOW_GLUT::
Request_Move(const int x,const int y)
{
    glutPositionWindow(x,y);
}
//#####################################################################
// Width
//#####################################################################
int OPENGL_WINDOW_GLUT::
Width() const
{
    return width;
}
//#####################################################################
// Height
//#####################################################################
int OPENGL_WINDOW_GLUT::
Height() const
{
    return height;
}
//#####################################################################
// Handle_Display_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Display_Glut()
{
    single->opengl_world.Render_World(false); // render, no selection
    glutSwapBuffers();
}
//#####################################################################
// Handle_Reshape_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Reshape_Glut(int w,int h)
{
    single->width=w;single->height=h;
    single->opengl_world.Handle_Reshape_Main();
}

//#####################################################################
// Handle_Special_Keypress_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Special_Keypress_Glut(int key,int x,int y)
{
    single->opengl_world.Handle_Keypress_Main(OPENGL_KEY::From_Glut_Special_Key(key,(glutGetModifiers()&GLUT_ACTIVE_CTRL)!=0,(glutGetModifiers()&GLUT_ACTIVE_ALT)!=0),x,y);
}
//#####################################################################
// Handle_Keypress_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Keypress_Glut(unsigned char key,int x,int y)
{
    if(single->opengl_world.prompt_mode) single->opengl_world.Handle_Keypress_Prompt(key);
    else single->opengl_world.Handle_Keypress_Main(OPENGL_KEY::From_Glut_Key(key,(glutGetModifiers()&GLUT_ACTIVE_CTRL)!=0,(glutGetModifiers()&GLUT_ACTIVE_ALT)!=0),x,y);
}
//#####################################################################
// Handle_Click_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Click_Glut(int button,int state,int x,int y)
{
    bool ctrl_pressed=(glutGetModifiers() & GLUT_ACTIVE_CTRL)!=0;
    bool shift_pressed=glutGetModifiers() & GLUT_ACTIVE_SHIFT;
    single->opengl_world.Handle_Click_Main(button,state,x,y,ctrl_pressed,shift_pressed);
}
//#####################################################################
// Handle_Drag_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Drag_Glut(int x,int y)
{
    single->opengl_world.Handle_Drag_Main(x,y);
}
//#####################################################################
// Handle_Idle_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Idle_Glut()
{
    single->opengl_world.Handle_Idle();
}
//#####################################################################
// Handle_Timer_Glut
//#####################################################################
void OPENGL_WINDOW_GLUT::
Handle_Timer_Glut(int value)
{
    single->opengl_world.Handle_Timer();
}
//#####################################################################
