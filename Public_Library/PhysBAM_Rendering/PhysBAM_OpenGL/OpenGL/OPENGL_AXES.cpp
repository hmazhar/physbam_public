//#####################################################################
// Copyright 2002-2005, Robert Bridson, Eran Guendelman, Eilene Hao, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_AXES
//##################################################################### 
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_AXES<T>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    glDisable(GL_LIGHTING);

    OpenGL_Begin(GL_LINES);

    if(in_color) glColor3f(1,.25,.25);
    else glColor3f(1,1,1);
    // draw all lines parallel to x axis
    if(draw_box){
        glVertex3f(box.min_corner.x,box.min_corner.y,box.min_corner.z);glVertex3f(box.max_corner.x,box.min_corner.y,box.min_corner.z);glVertex3f(box.min_corner.x,box.max_corner.y,box.min_corner.z);glVertex3f(box.max_corner.x,box.max_corner.y,box.min_corner.z);
        glVertex3f(box.min_corner.x,box.max_corner.y,box.max_corner.z);glVertex3f(box.max_corner.x,box.max_corner.y,box.max_corner.z);glVertex3f(box.min_corner.x,box.min_corner.y,box.max_corner.z);glVertex3f(box.max_corner.x,box.min_corner.y,box.max_corner.z);}
    else{glVertex3f(box.min_corner.x,0,0);glVertex3f(box.max_corner.x,0,0);}
    if(draw_xz_grid) for(T z=box.min_corner.z;z<=box.max_corner.z;z+=grid_spacing){glVertex3f(box.min_corner.x,box.min_corner.y,z);glVertex3f(box.max_corner.x,box.min_corner.y,z);}
    if(draw_xy_grid) for(T y=box.min_corner.y;y<=box.max_corner.y;y+=grid_spacing){glVertex3f(box.min_corner.x,y,box.min_corner.z);glVertex3f(box.max_corner.x,y,box.min_corner.z);}

    if(in_color) glColor3f(.25f,1,.25f);
    // draw all lines parallel to y axis
    if(draw_box){
        glVertex3f(box.min_corner.x,box.min_corner.y,box.min_corner.z);glVertex3f(box.min_corner.x,box.max_corner.y,box.min_corner.z);glVertex3f(box.max_corner.x,box.min_corner.y,box.min_corner.z);glVertex3f(box.max_corner.x,box.max_corner.y,box.min_corner.z);
        glVertex3f(box.max_corner.x,box.min_corner.y,box.max_corner.z);glVertex3f(box.max_corner.x,box.max_corner.y,box.max_corner.z);glVertex3f(box.min_corner.x,box.min_corner.y,box.max_corner.z);glVertex3f(box.min_corner.x,box.max_corner.y,box.max_corner.z);}
    else{
        glVertex3f(0,box.min_corner.y,0);glVertex3f(0,box.max_corner.y,0);}
    if(draw_yz_grid) for(T z=box.min_corner.z;z<=box.max_corner.z;z+=grid_spacing){glVertex3f(box.min_corner.x,box.min_corner.y,z);glVertex3f(box.min_corner.x,box.max_corner.y,z);}
    if(draw_xy_grid) for(T x=box.min_corner.x;x<=box.max_corner.x;x+=grid_spacing){glVertex3f(x,box.min_corner.y,box.min_corner.z);glVertex3f(x,box.max_corner.y,box.min_corner.z);}

    if(in_color) glColor3f(.25f,.25f,1);
    // draw all lines parallel to z axis
    if(draw_box){
        glVertex3f(box.min_corner.x,box.min_corner.y,box.min_corner.z);glVertex3f(box.min_corner.x,box.min_corner.y,box.max_corner.z);glVertex3f(box.max_corner.x,box.min_corner.y,box.min_corner.z);glVertex3f(box.max_corner.x,box.min_corner.y,box.max_corner.z);
        glVertex3f(box.max_corner.x,box.max_corner.y,box.min_corner.z);glVertex3f(box.max_corner.x,box.max_corner.y,box.max_corner.z);glVertex3f(box.min_corner.x,box.max_corner.y,box.min_corner.z);glVertex3f(box.min_corner.x,box.max_corner.y,box.max_corner.z);}
    else{glVertex3f(0,0,box.min_corner.z);glVertex3f(0,0,box.max_corner.z);}
    if(draw_yz_grid) for(T y=box.min_corner.y;y<=box.max_corner.y;y+=grid_spacing){glVertex3f(box.min_corner.x,y,box.min_corner.z);glVertex3f(box.min_corner.x,y,box.max_corner.z);}
    if(draw_xz_grid) for(T x=box.min_corner.x;x<=box.max_corner.x;x+=grid_spacing){glVertex3f(x,box.min_corner.y,box.min_corner.z);glVertex3f(x,box.min_corner.y,box.max_corner.z);}

    OpenGL_End();
    glPopMatrix();
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_AXES<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(box));
}
//#####################################################################
template class OPENGL_AXES<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_AXES<double>;
#endif
