//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_DEBUG_PARTICLES_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_DEBUG_PARTICLES_2D<T>::
OPENGL_DEBUG_PARTICLES_2D(GEOMETRY_PARTICLES<TV>& particle_input,const OPENGL_COLOR& color_input)
    :particles(particle_input),default_color(color_input),velocity_color(OPENGL_COLOR(1,(T).078,(T).576)),draw_velocities(false),draw_arrows(true),scale_velocities((T).025)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_DEBUG_PARTICLES_2D<T>::
~OPENGL_DEBUG_PARTICLES_2D()
{
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_DEBUG_PARTICLES_2D<T>::
Use_Bounding_Box() const
{
    return particles.X.Size()>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_DEBUG_PARTICLES_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>::Bounding_Box(particles.X));
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
    glPointSize(5);
    glDisable(GL_LIGHTING);

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    ARRAY_VIEW<VECTOR<T,3> >* colors=particles.array_collection->template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    ARRAY_VIEW<T>* sizes=particles.array_collection->template Get_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);
    ARRAY_VIEW<TV>* V=particles.array_collection->template Get_Array<TV>(ATTRIBUTE_ID_V);

    if(draw_velocities && V && mode!=GL_SELECT){
        glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int i=1;i<=particles.X.m;i++){
            TV X=particles.X(i);
            TV Y=X+(*V)(i)*scale_velocities;
            if(draw_arrows) OPENGL_SHAPES::Draw_Arrow(X,Y);
            else OpenGL_Line(X,Y);}
        OpenGL_End();
        glPopAttrib();}

    if(mode==GL_SELECT) glPushName(0);
    for(int i=1;i<=particles.X.m;i++){
        if(mode==GL_SELECT) glLoadName(i);

        if(colors) OPENGL_COLOR((*colors)(i)).Send_To_GL_Pipeline();
        else default_color.Send_To_GL_Pipeline();

        if(sizes) OPENGL_SHAPES::Draw_Circle(particles.X(i),(*sizes)(i),20,true);
        else{
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(particles.X(i));
            OpenGL_End();}}
    if(mode==GL_SELECT) glPopName();

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Select_Point
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Select_Point(int index)
{
    //Set_Point_Color(index,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Select_Points
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Select_Points(const ARRAY<int> &indices)
{
    //Set_Point_Colors(indices,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Clear_Selection
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Clear_Selection()
{
    //Store_Point_Colors(false);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_DEBUG_PARTICLES_2D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    LOG::cout<<"Get_Selection "<<buffer_size<<std::endl;
    if(buffer_size==1){
        OPENGL_SELECTION_DEBUG_PARTICLES_2D<T> *selection=new OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>(this);
        selection->index=buffer[0];
        return selection;}
    else return 0;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type!=OPENGL_SELECTION::DEBUG_PARTICLES_2D) return;
    OPENGL_SELECTION_DEBUG_PARTICLES_2D<T> *real_selection=(OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>*)selection;
    Select_Point(real_selection->index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Clear_Highlight()
{
    Clear_Selection();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    return object->Bounding_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    if(selection->type!=OPENGL_SELECTION::DEBUG_PARTICLES_2D) return;
    output_stream<<"Particle "<<dynamic_cast<OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>&>(*selection).index<<std::endl;
}

template class OPENGL_DEBUG_PARTICLES_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_DEBUG_PARTICLES_2D<double>;
#endif
