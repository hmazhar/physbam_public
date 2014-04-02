//#####################################################################
// Copyright 2004-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SYMMETRIC_MATRIX_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT|GL_TEXTURE_BIT|GL_LINE_BIT);
    glLineWidth(1);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
    OpenGL_Begin(GL_LINES);
    for(int i=1;i<=entries.m;i++){
        VECTOR<T,3> node=entries(i).x;MATRIX<T,3> line=size*entries(i).y;VECTOR<bool,3> p=entries(i).z;
        (p.x?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(1),node+line.Column(1));
        (p.y?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(2),node+line.Column(2));
        (p.z?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(3),node+line.Column(3));}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T>::
Update()
{
    if(!field.counts.x || !field.counts.y || !field.counts.z){entries.Clean_Memory();return;}

    int m_start=1,m_end=grid.counts.x,n_start=1,n_end=grid.counts.y,mn_start=1,mn_end=grid.counts.z;
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    if(slice && slice->mode!=OPENGL_SLICE::NO_SLICE){
        VECTOR<int,3> domain_start(m_start,n_start,mn_start),domain_end(m_end,n_end,mn_end);
        if((slice->mode == OPENGL_SLICE::CELL_SLICE && (grid.MAC_offset==0 || slice->index < domain_start[slice->axis] || slice->index > domain_end[slice->axis])) ||
           (slice->mode == OPENGL_SLICE::NODE_SLICE && (grid.MAC_offset==0.5 || slice->index < domain_start[slice->axis] || slice->index > domain_end[slice->axis]))) return;
        switch(slice->axis){
            case 1:m_start=m_end=slice->index;break;
            case 2:n_start=n_end=slice->index;break;
            case 3:mn_start=mn_end=slice->index;break;}}

    int count=0;
    for(int i=m_start;i<=m_end;i++)for(int j=m_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) if(field(i,j,ij).Frobenius_Norm_Squared())count++;
    entries.Resize(count);count=0;
    DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> U;
    for(int i=m_start;i<=m_end;i++)for(int j=m_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++) if(field(i,j,ij).Frobenius_Norm_Squared()){
        field(i,j,ij).Solve_Eigenproblem(D,U);entries(++count)=TRIPLE<VECTOR<T,3>,MATRIX<T,3>,VECTOR<bool,3> >(grid.X(i,j,ij),U*D,VECTOR<bool,3>(D.x11>0,D.x22>0,D.x33>0));}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.domain);
}
//#####################################################################
template class OPENGL_SYMMETRIC_MATRIX_FIELD_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SYMMETRIC_MATRIX_FIELD_3D<double>;
#endif
