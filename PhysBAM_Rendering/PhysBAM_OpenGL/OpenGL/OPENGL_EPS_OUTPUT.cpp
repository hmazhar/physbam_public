//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_EPS_OUTPUT.h>
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
using namespace PhysBAM;
template<class T> OPENGL_EPS_OUTPUT<T>::
OPENGL_EPS_OUTPUT(const std::string& filename)
    :stream(*FILE_UTILITIES::Safe_Open_Output(filename,false,false))
{
    Head();
}
template<class T> OPENGL_EPS_OUTPUT<T>::
~OPENGL_EPS_OUTPUT()
{
    Tail();
    delete &stream;
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Begin(int mode)
{
    glmode=mode;
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
End()
{
    Transform_Buffer();
    switch(glmode){
        case GL_POINTS: for(int i=1;i<=buffer.m;i++) Draw_Point(buffer(i)); break;
        case GL_LINES: for(int i=1;i<=buffer.m;i+=2) Draw_Line(buffer(i),buffer(i+1)); break;
        case GL_TRIANGLES: for(int i=1;i<=buffer.m;i+=3) Draw_Polygon(i,3); break;
        case GL_QUADS: for(int i=1;i<=buffer.m;i+=4) Draw_Polygon(i,4); break;
        case GL_POLYGON: Draw_Polygon(1,buffer.m); break;
        case GL_LINE_LOOP: for(int i=2;i<=buffer.m;i++) Draw_Line(buffer.Last(),buffer(1)); break;
        default:LOG::cout<<"Unhandled mode "<<glmode<<std::endl;break;
    }
    buffer.Remove_All();
}
template<class T> template<class T2> void OPENGL_EPS_OUTPUT<T>::
Vertex(const VECTOR<T2,3>& p)
{
    buffer.Append(VECTOR<T,3>(p));
}
template<class T> template<class T2> void OPENGL_EPS_OUTPUT<T>::
Vertex(const VECTOR<T2,2>& p)
{
    buffer.Append(VECTOR<T,3>(VECTOR<T,2>(p)));
}
template<class T> template<class T2> void OPENGL_EPS_OUTPUT<T>::
Vertex(const VECTOR<T2,1>& p)
{
    buffer.Append(VECTOR<T,3>(p.x,0,0));
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Head()
{
    float view[4];
    glGetFloatv(GL_VIEWPORT,view);

    Emit("%!PS-Adobe-3.0 EPSF-3.0\n");
    Emit("%%BoundingBox: ");
    stream<<view[0]<<" "<<view[1]<<" "<<(view[2]+view[0])<<" "<<(view[3]+view[1])<<"\n";
    Emit("0 setlinewidth\n");
    Emit("/pointradius .1 def\n");
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Tail()
{
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Emit(const char* s)
{
    stream<<s;
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Emit(const TV& p)
{
    stream<<p.x<<" "<<p.y<<" ";
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Draw_Point(const TV& p)
{
    Emit(p);
    Emit("pointradius 0 360 arc stroke\n");
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Draw_Line(const TV& a,const TV& b)
{
    Emit(a);
    Emit("moveto ");
    Emit(b);
    Emit("lineto stroke\n");
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Draw_Polygon(int i,int n)
{
    GLint mode[2]={0};
    glGetIntegerv(GL_POLYGON_MODE,mode);
    const char* str = "closepath stroke\n";
    if(mode[0] == GL_FILL) str = "closepath fill stroke\n";
    Emit(buffer(i));
    Emit("moveto ");
    for(int j=i+1;j<i+n;j++){
        Emit(buffer(j));
        Emit("lineto ");}
    Emit(str);
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Set_Color(const TV& color)
{
    stream<<color.x<<" "<<color.y<<" "<<color.z<<" setrgbcolor"<<std::endl;
}
template<class T> void OPENGL_EPS_OUTPUT<T>::
Transform_Buffer()
{
    MATRIX<T,4> proj,model;
    T view[4];

    if(sizeof(T)==sizeof(float)){
        glGetFloatv(GL_PROJECTION_MATRIX,(GLfloat*)proj.x);
        glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)model.x);
        glGetFloatv(GL_VIEWPORT,(GLfloat*)view);}
    else{
        glGetDoublev(GL_PROJECTION_MATRIX,(GLdouble*)proj.x);
        glGetDoublev(GL_MODELVIEW_MATRIX,(GLdouble*)model.x);
        glGetDoublev(GL_VIEWPORT,(GLdouble*)view);}

    for(int i=1;i<=buffer.m;i++){
        VECTOR<T,4> obj=buffer(i).Append(1);
        VECTOR<T,4> eye=model*obj;
        VECTOR<T,4> clip=proj*eye;
        TV device=clip.Remove_Index(4)/clip(4);
        TV window((device.x+1)*(view[2]/2)+view[0],(device.y+1)*(view[3]/2)+view[1],device.z);
        buffer(i)=window;}
}
template class OPENGL_EPS_OUTPUT<float>;
template void OPENGL_EPS_OUTPUT<float>::Vertex<float>(VECTOR<float,2> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<float>(VECTOR<float,3> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<double>(VECTOR<double,2> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<double>(VECTOR<double,3> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_EPS_OUTPUT<double>;
template void OPENGL_EPS_OUTPUT<double>::Vertex<float>(VECTOR<float,2> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<float>(VECTOR<float,3> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<double>(VECTOR<double,2> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<double>(VECTOR<double,3> const&);
#endif
