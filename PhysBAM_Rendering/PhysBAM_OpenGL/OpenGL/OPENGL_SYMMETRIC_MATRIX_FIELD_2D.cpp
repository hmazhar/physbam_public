//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SYMMETRIC_MATRIX_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SYMMETRIC_MATRIX_FIELD_2D<T>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT|GL_TEXTURE_BIT|GL_LINE_BIT);
    glLineWidth(1);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
    OpenGL_Begin(GL_LINES);
    for(int i=lines.domain.min_corner.x;i<=lines.domain.max_corner.x;i++)for(int j=lines.domain.min_corner.y;j<=lines.domain.max_corner.y;j++){
        VECTOR<T,2> node=grid.Node(i,j);MATRIX<T,2> line=size*lines(i,j);
        (positive(i,j).x?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(1),node+line.Column(1));
        (positive(i,j).y?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(2),node+line.Column(2));}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_SYMMETRIC_MATRIX_FIELD_2D<T>::
Update()
{
    lines.Resize(grid.Domain_Indices(1-field.domain.min_corner.x));positive.Resize(grid.Domain_Indices(1-field.domain.min_corner.x));
    DIAGONAL_MATRIX<T,2> D;MATRIX<T,2> U;
    for(int i=lines.domain.min_corner.x;i<=lines.domain.max_corner.x;i++)for(int j=lines.domain.min_corner.y;j<=lines.domain.max_corner.y;j++){
        field(i,j).Solve_Eigenproblem(D,U);lines(i,j)=U*D;positive(i,j)=PAIR<bool,bool>(D.x11>0,D.x22>0);}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SYMMETRIC_MATRIX_FIELD_2D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y,grid.domain.max_corner.y,0,0);
}
//#####################################################################
template class OPENGL_SYMMETRIC_MATRIX_FIELD_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SYMMETRIC_MATRIX_FIELD_2D<double>;
#endif
