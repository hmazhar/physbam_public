//#####################################################################
// Copyright 2002,2003,Neil Molino,Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// File OPENGL_SHAPES.cpp
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
using namespace PhysBAM::OPENGL_SHAPES;

template<class T>
void OPENGL_SHAPES::Draw_Vector(const VECTOR<T,2>& from,const VECTOR<T,2>& v,OPENGL_COLOR color,const float size)
{
    if(v.Magnitude()==0){Draw_Dot(from,color,size);return;}
    glMatrixMode(GL_MODELVIEW); glPushMatrix();
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT | GL_POINT_BIT);
    glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);glLineWidth(size);glPointSize(7); 
    color.Send_To_GL_Pipeline();
    
    T arrow_relative_axial_length=.1,arrow_axial_length=arrow_relative_axial_length*v.Magnitude();
    VECTOR<T,2> to=from+v,a=from+(1-arrow_relative_axial_length)*v;
    
    VECTOR<T,2> lateral_direction=v.Rotate_Clockwise_90();
    lateral_direction.Normalize();
    T arrow_lateral_length=arrow_axial_length*.5;
    VECTOR<T,2> a1=a+arrow_lateral_length*lateral_direction,
                                    a2=a-arrow_lateral_length*lateral_direction;
    
    OpenGL_Begin(GL_LINES);
    OpenGL_Line(from,to);
    OpenGL_Line(a1,to);
    OpenGL_Line(a2,to);
    OpenGL_End();
    glPopAttrib();
    glPopMatrix();
}

template<class T>
void OPENGL_SHAPES::Draw_Vector(const VECTOR<T,3>& from,const VECTOR<T,3>& v,OPENGL_COLOR color,const float size)
{
    if(v.Magnitude()==0){Draw_Dot(from,color,size);return;}
    glMatrixMode(GL_MODELVIEW); glPushMatrix();
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT | GL_POINT_BIT);
    glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);glLineWidth(size);glPointSize(7); 
    color.Send_To_GL_Pipeline();
    
    T arrow_relative_axial_length=.1,arrow_axial_length=arrow_relative_axial_length*v.Magnitude();
    VECTOR<T,3> to=from+v,a=from+(1-arrow_relative_axial_length)*v;
    
    VECTOR<T,3> lateral_direction = v.Orthogonal_Vector();
    /*VECTOR<T,3>::Cross_Product(VECTOR<T,3>(1,1,1),v); 
    if(lateral_direction.Magnitude()==0) lateral_direction = VECTOR<T,3>::Cross_Product(VECTOR<T,3>(1,0,0),v); */
    lateral_direction.Normalize();
    T arrow_lateral_length=arrow_axial_length*.5;
    VECTOR<T,3> a1=a+arrow_lateral_length*lateral_direction,
                                    a2=a-arrow_lateral_length*lateral_direction;
    
    OpenGL_Begin(GL_LINES);
    OpenGL_Line(from,to);
    OpenGL_Line(a1,to);
    OpenGL_Line(a2,to);
    OpenGL_End();
    glPopAttrib();
    glPopMatrix();
}

template void OPENGL_SHAPES::Draw_Vector<float>(const VECTOR<float,2>&,const VECTOR<float,2>&,OPENGL_COLOR,const float);
template void OPENGL_SHAPES::Draw_Vector<float>(const VECTOR<float,3>&,const VECTOR<float,3>&,OPENGL_COLOR,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void OPENGL_SHAPES::Draw_Vector<double>(const VECTOR<double,2>&,const VECTOR<double,2>&,OPENGL_COLOR,const float);
template void OPENGL_SHAPES::Draw_Vector<double>(const VECTOR<double,3>&,const VECTOR<double,3>&,OPENGL_COLOR,const float);
#endif
