//#####################################################################
// Copyright 2004, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    if(draw_cells) OPENGL_SCALAR_FIELD_2D<T,T>::Display(in_color);
    if(draw_normals){
        levelset.Compute_Normals();
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);OPENGL_COLOR::White().Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int i=1;i<=levelset.grid.counts.x;i++) for(int j=1;j<=levelset.grid.counts.y;j++) if(!active_cells || (*active_cells)(i,j)){
            OpenGL_Line(grid.X(i,j),grid.X(i,j)+T(0.01)*(*levelset.normals)(i,j));}
        OpenGL_End();
        glPopAttrib();}
    if(draw_area&&opengl_triangulated_area){glDepthMask(GL_FALSE);opengl_triangulated_area->Display(in_color);glDepthMask(GL_TRUE);}
    if(draw_curve&&opengl_segmented_curve_2d){opengl_segmented_curve_2d->Display(in_color);}
    glPopMatrix();
}
//#####################################################################
// Function Set_Inside_And_Outside_Colors
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Set_Inside_And_Outside_Colors(const OPENGL_COLOR& inside_color_input,const OPENGL_COLOR& outside_color_input)
{
    inside_color=inside_color_input;outside_color=outside_color_input;
    this->color_maps(1);
    this->color_maps(1)=OPENGL_COLOR_RAMP<T>::Levelset_Color_Constant_Ramp(inside_color,outside_color);
}
//#####################################################################
// Function Toggle_Normals
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Toggle_Normals()
{
    if(draw_normals)draw_normals=false;else draw_normals=true;
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_LEVELSET_2D<T>::
Update()
{
    OPENGL_SCALAR_FIELD_2D<T,T>::Update();
    if(levelset.phi.counts.x>1&&levelset.phi.counts.y>1){
        if(opengl_triangulated_area){
            delete &(opengl_triangulated_area->triangulated_area);delete opengl_triangulated_area;opengl_triangulated_area=0;}
        if(opengl_segmented_curve_2d){
            delete &(opengl_segmented_curve_2d->curve);delete opengl_segmented_curve_2d;opengl_segmented_curve_2d=0;}
        if(draw_area||draw_curve){
            DUALCONTOUR_2D<T> dualcontour(levelset);dualcontour.Dualcontour();
            OPENGL_COLOR color=(dominant_sign==1)?OPENGL_COLOR::Red():OPENGL_COLOR::Blue();
            opengl_triangulated_area=new OPENGL_TRIANGULATED_AREA<T>(*dualcontour.Get_Triangulated_Area(dominant_sign),false,OPENGL_COLOR::Red(),OPENGL_COLOR::Blue(),color);
            opengl_segmented_curve_2d=new OPENGL_SEGMENTED_CURVE_2D<T>(*dualcontour.Get_Segmented_Curve());}}
}
//#####################################################################
template class OPENGL_LEVELSET_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_LEVELSET_2D<double>;
#endif
