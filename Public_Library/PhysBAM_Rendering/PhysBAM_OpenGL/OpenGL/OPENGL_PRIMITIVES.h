//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Robert Bridson, Sergey Koltakov, Neil Molino, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace OPENGL_PRIMITIVES
//#####################################################################
#ifndef __OPENGL_PRIMITIVES__
#define __OPENGL_PRIMITIVES__
#include <cstdio>
#ifdef _WIN32
#include <windows.h>
#endif

#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_EPS_OUTPUT.h>
namespace PhysBAM{

#define WANT_OPENGL_EPS_OUTPUT // Uncomment this if you want support for dumping eps.
extern OPENGL_EPS_OUTPUT<float>* opengl_eps_output;
#ifdef WANT_OPENGL_EPS_OUTPUT
#define IF_OPENGL_EPS_OUTPUT(x) if(opengl_eps_output){x;}
#else
#define IF_OPENGL_EPS_OUTPUT(x)
#endif

inline void OpenGL_Begin(GLenum mode)
{glBegin(mode);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));}

inline void OpenGL_End()
{glEnd();IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Eps_Emit(const char* str)
{IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Emit(str));}

inline void OpenGL_Rotate(const ROTATION<VECTOR<float,3> >& r)
{float angle;VECTOR<float,3>axis;r.Get_Angle_Axis(angle,axis);glRotatef(angle*180/(float)pi,axis.x,axis.y,axis.z);}

inline void OpenGL_Rotate(const ROTATION<VECTOR<double,3> >& r)
{double angle;VECTOR<double,3> axis;r.Get_Angle_Axis(angle,axis);glRotated(angle*180/pi,axis.x,axis.y,axis.z);}

inline void OpenGL_Translate(const VECTOR<float,2>& v)
{glTranslatef(v.x,v.y,(float)0);}

inline void OpenGL_Translate(const VECTOR<double,2>& v)
{glTranslated(v.x,v.y,(double)0);}

inline void OpenGL_Translate(const VECTOR<float,3>& v)
{glTranslatef(v.x,v.y,v.z);}

inline void OpenGL_Translate(const VECTOR<double,3>& v)
{glTranslated(v.x,v.y,v.z);}

inline void OpenGL_Scale(const VECTOR<float,3>& v)
{glScalef(v.x,v.y,v.z);}

inline void OpenGL_Scale(const VECTOR<double,3>& v)
{glScaled(v.x,v.y,v.z);}

template<class T>
inline void OpenGL_Transform(const FRAME<VECTOR<T,3> >& frame)
{OpenGL_Translate(frame.t);OpenGL_Rotate(frame.r.Normalized());}

inline void OpenGL_Vertex(const VECTOR<float,3>& v)
{glVertex3f(v.x,v.y,v.z);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<float,2>& v)
{glVertex2f(v.x,v.y);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<float,1>& v)
{glVertex2f(v.x,(float)0);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<float,2>(v.x,0)));}

inline void OpenGL_Normal(const VECTOR<float,3>& n)
{glNormal3f(n.x,n.y,n.z);}

inline void OpenGL_RasterPos(const VECTOR<float,1>& v)
{glRasterPos2f(v.x,(float)0);}

inline void OpenGL_RasterPos(const VECTOR<float,2>& v)
{glRasterPos2f(v.x,v.y);}

inline void OpenGL_RasterPos(const VECTOR<float,3>& v)
{glRasterPos3f(v.x,v.y,v.z);}

inline void OpenGL_Vertex(const VECTOR<double,3>& v)
{glVertex3d(v.x,v.y,v.z);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<double,2>& v)
{glVertex2d(v.x,v.y);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<double,1>& v)
{glVertex2d(v.x,(double)0);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<double,2>(v.x,0)));}

template<class T,int d>
inline void OpenGL_Line(const VECTOR<T,d>& a, const VECTOR<T,d>& b)
{OpenGL_Vertex(a);OpenGL_Vertex(b);}

template<class T,int d>
inline void OpenGL_Triangle(const VECTOR<T,d>& a, const VECTOR<T,d>& b, const VECTOR<T,d>& c)
{OpenGL_Vertex(a);OpenGL_Vertex(b);OpenGL_Vertex(c);}

inline void OpenGL_Normal(const VECTOR<double,3>& n)
{glNormal3d(n.x,n.y,n.z);}

inline void OpenGL_RasterPos(const VECTOR<double,1>& v)
{glRasterPos2f(v.x,(double)0);}

inline void OpenGL_RasterPos(const VECTOR<double,2>& v)
{glRasterPos2d(v.x,v.y);}

inline void OpenGL_RasterPos(const VECTOR<double,3>& v)
{glRasterPos3d(v.x,v.y,v.z);}

inline void OpenGL_LookFrom(const FRAME<VECTOR<double,3> >& frame)
{OpenGL_Rotate(frame.r.Inverse().Normalized());OpenGL_Translate(-frame.t);}

template<class TV>
inline void OpenGL_String(const TV& position,const std::string& str,void* font=GLUT_BITMAP_HELVETICA_12)
{OpenGL_RasterPos(position);
for(unsigned int j=0;j<str.length();j++) glutBitmapCharacter(font,str[j]);}

inline void OpenGL_Quad_2D(const VECTOR<float,2> &bottom_left,const VECTOR<float,2> &top_right)
{glVertex2f(bottom_left.x,bottom_left.y);glVertex2f(bottom_left.x,top_right.y);
 glVertex2f(top_right.x,top_right.y);glVertex2f(top_right.x,bottom_left.y);
 IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<float,2>(bottom_left.x,bottom_left.y));opengl_eps_output->Vertex(VECTOR<float,2>(bottom_left.x,top_right.y));
     opengl_eps_output->Vertex(VECTOR<float,2>(top_right.x,top_right.y));opengl_eps_output->Vertex(VECTOR<float,2>(top_right.x,bottom_left.y)));}

inline void OpenGL_Quad_2D(const VECTOR<double,2> &bottom_left,const VECTOR<double,2> &top_right)
{glVertex2d(bottom_left.x,bottom_left.y);glVertex2d(bottom_left.x,top_right.y);
 glVertex2d(top_right.x,top_right.y);glVertex2d(top_right.x,bottom_left.y);
 IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<double,2>(bottom_left.x,bottom_left.y));opengl_eps_output->Vertex(VECTOR<double,2>(bottom_left.x,top_right.y));
     opengl_eps_output->Vertex(VECTOR<double,2>(top_right.x,top_right.y));opengl_eps_output->Vertex(VECTOR<double,2>(top_right.x,bottom_left.y)));}

template<class T>
inline void OpenGL_Quad(const VECTOR<T,3> &bottom_left,const VECTOR<T,3> &right,const VECTOR<T,3>& up)
{OpenGL_Vertex(bottom_left);OpenGL_Vertex(bottom_left+up);OpenGL_Vertex(bottom_left+up+right);OpenGL_Vertex(bottom_left+right);}

template<class T>
inline void OpenGL_Clip_Plane(GLenum id,const PLANE<T> &plane)
{GLdouble equation[4]={plane.normal.x,plane.normal.y,plane.normal.z,-VECTOR<T,3>::Dot_Product(plane.normal,plane.x1)};
glClipPlane(id,equation);}

}
#endif
