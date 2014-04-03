//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
namespace PhysBAM{
//#####################################################################
// OPENGL_FACE_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
OPENGL_FACE_SCALAR_FIELD_1D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<1> > &face_values_input,OPENGL_COLOR point_color_input,OPENGL_COLOR line_color_input)
:grid(grid_input),face_values(face_values_input),x_face_values(face_values.Component(1)),point_color(point_color_input),line_color(line_color_input),scale(1)
{
}
//#####################################################################
// ~OPENGL_FACE_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
~OPENGL_FACE_SCALAR_FIELD_1D()
{
}
//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    line_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_STRIP);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        OpenGL_Vertex(VECTOR<T,2>(iterator.Location().x,scale*x_face_values(iterator.Face_Index())));}
    OpenGL_End();
    glColor3f(0,1,1);
    glPointSize(3.0);
    point_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_POINTS);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        OpenGL_Vertex(VECTOR<T,2>(iterator.Location().x,scale*x_face_values(iterator.Face_Index())));}
    OpenGL_End();
    glPopAttrib();
}
template<class T> void
Display_Bool_Helper(const OPENGL_FACE_SCALAR_FIELD_1D<T,bool>& self,const int in_color)
{
    typedef typename GRID<VECTOR<T,1> >::FACE_ITERATOR FACE_ITERATOR;
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glPointSize(8.0);
    self.point_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_POINTS);
    for(FACE_ITERATOR iterator(self.grid);iterator.Valid();iterator.Next()) if(self.x_face_values(iterator.Face_Index())){
        OpenGL_Vertex(VECTOR<T,2>(iterator.Location().x,(T)0));}
    OpenGL_End();
    glPopAttrib();
}
template<> void OPENGL_FACE_SCALAR_FIELD_1D<float,bool>::
Display(const int in_color) const
{
    Display_Bool_Helper(*this,in_color);
}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> void OPENGL_FACE_SCALAR_FIELD_1D<double,bool>::
Display(const int in_color) const
{
    Display_Bool_Helper(*this,in_color);
}
#endif
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(grid.domain.min_corner.x,grid.domain.max_corner.x,0,0,0,0));
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    // TODO: this should also interpolate to particles
    if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_1D && grid.Is_MAC_Grid()){
        VECTOR<int,1> index=((OPENGL_SELECTION_GRID_CELL_1D<T>*)selection)->index;
        T2 left=x_face_values(index.x);
        T2 right=x_face_values(index.x+1);
        output_stream<<"    left = "<<left<<",right = "<<right<<std::endl;}
}
//#####################################################################
template class OPENGL_FACE_SCALAR_FIELD_1D<float,int>;
template class OPENGL_FACE_SCALAR_FIELD_1D<float,bool>;
template class OPENGL_FACE_SCALAR_FIELD_1D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_FACE_SCALAR_FIELD_1D<double,int>;
template class OPENGL_FACE_SCALAR_FIELD_1D<double,bool>;
template class OPENGL_FACE_SCALAR_FIELD_1D<double,double>;
#endif
}
