//#####################################################################
// Copyright 2006-2009, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_SCALAR_FIELD_1D<T,T2>::
OPENGL_SCALAR_FIELD_1D(const GRID<TV>& grid,ARRAY<T2,VECTOR<int,1> >& values,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :grid(grid),values(values),point_color(point_color),line_color(line_color),scale(1)
{}
//#####################################################################
// Function OPENGL_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_SCALAR_FIELD_1D<T,T2>::
~OPENGL_SCALAR_FIELD_1D()
{}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_SCALAR_FIELD_1D<T,T2>::
Display(const int in_color) const
{
    int min_corner=values.Domain_Indices().min_corner.x;
    int max_corner=values.Domain_Indices().max_corner.x;

    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    line_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_STRIP);
    for(int i=min_corner;i<=max_corner;i++) OpenGL_Vertex(VECTOR<T,3>(grid.Axis_X(i,1),values(i)*scale,(T)0));
    OpenGL_End();
    glColor3f(0,1,1);
    glPointSize(3.0);
    point_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_POINTS);
    for(int i=min_corner;i<=max_corner;i++) OpenGL_Vertex(VECTOR<T,3>(grid.Axis_X(i,1),values(i)*scale,(T)0));
    OpenGL_End();
    glPopAttrib();
}
template<class T>
void Display_Bool_Helper(const OPENGL_SCALAR_FIELD_1D<T,bool>& self,const int in_color)
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(0,1,1);
    glPointSize(8.0);
    self.point_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_POINTS);
    for(int i=self.values.domain.min_corner.x;i<=self.values.domain.max_corner.x;i++) if(self.values(i)) OpenGL_Vertex(VECTOR<T,3>(self.grid.Axis_X(i,1),(T)0,(T)0));
    OpenGL_End();
    glPopAttrib();
}
namespace PhysBAM{
template<> void OPENGL_SCALAR_FIELD_1D<float,bool>::
Display(const int in_color) const
{
    Display_Bool_Helper(*this,in_color);
}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> void OPENGL_SCALAR_FIELD_1D<double,bool>::
Display(const int in_color) const
{
    Display_Bool_Helper(*this,in_color);
}
#endif
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_SCALAR_FIELD_1D<T,T2>::
Bounding_Box() const
{
    int min_corner=values.Domain_Indices().min_corner.x;
    int max_corner=values.Domain_Indices().max_corner.x;
    RANGE<VECTOR<float,3> > box=RANGE<VECTOR<float,3> >::Empty_Box();
    for(int i=min_corner;i<=max_corner;i++) box.Enlarge_To_Include_Point(VECTOR<float,3>(grid.Axis_X(i,1),values(i),(T)0));
    LOG::cout<<"box is "<<box<<std::endl;
    return box;
}
//#####################################################################
template class OPENGL_SCALAR_FIELD_1D<float,float>;
template class OPENGL_SCALAR_FIELD_1D<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SCALAR_FIELD_1D<double,double>;
template class OPENGL_SCALAR_FIELD_1D<double,bool>;
#endif
