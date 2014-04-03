//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Sergey Koltakov, Neil Molino, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_SEGMENTED_CURVE_2D<T>::
OPENGL_SEGMENTED_CURVE_2D(const SEGMENTED_CURVE_2D<T>& curve_input,const OPENGL_COLOR &color_input)
    :curve(curve_input),color(color_input),color_gray(color_input.Grayscale()),
    vertex_color(OPENGL_COLOR::Green(0.9)),vertex_position_color(OPENGL_COLOR::Magenta()),velocity_color(OPENGL_COLOR::Cyan()),
    draw_vertices(false),draw_vertex_positions(false),draw_velocities(false),velocity_scale(0.025),current_selection(0)
{}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    if(in_color) color.Send_To_GL_Pipeline();
    else color_gray.Send_To_GL_Pipeline();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);
    if (mode == GL_SELECT)
    {
        glPushName(1);
        Draw_Vertices_For_Selection();
        glLoadName(2);
        Draw_Segments_For_Selection();
        glPopName();
    }
    else
    {
        OpenGL_Begin(GL_LINES);
        for(int t=1; t<=curve.mesh.elements.m; t++){
            int i=curve.mesh.elements(t)(1),j=curve.mesh.elements(t)(2);
            OpenGL_Line(curve.particles.X(i),curve.particles.X(j));}
        OpenGL_End();

        if (current_selection) {
            if (current_selection->type == OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D) {
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)current_selection)->index;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(curve.particles.X(index));
            }
            else if (current_selection->type == OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D) {
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)current_selection)->index;
                int node1,node2;curve.mesh.elements(index).Get(node1,node2);
                OPENGL_SELECTION::Draw_Highlighted_Segment(curve.particles.X(node1),curve.particles.X(node2));
            }
        }
    }

    if (draw_vertices) {
        vertex_color.Send_To_GL_Pipeline();
        glPointSize(5.0f);
        OpenGL_Begin(GL_POINTS);
        for(int t=1; t<=curve.particles.array_collection->Size(); t++){
            OpenGL_Vertex(curve.particles.X(t));
        }
        OpenGL_End();
    }

    if (draw_vertex_positions) {
        vertex_position_color.Send_To_GL_Pipeline();
        for(int t=1; t<=curve.particles.array_collection->Size(); t++){
            VECTOR<float,3> world_space_pos=World_Space_Point(VECTOR<float,2>(curve.particles.X(t)));
            OpenGL_String(curve.particles.X(t),STRING_UTILITIES::string_sprintf("<%f %f>",world_space_pos.x,world_space_pos.y));
        }
    }

    if (draw_velocities && curve.particles.store_velocity){
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int t=1;t<=curve.particles.array_collection->Size();t++)
            OPENGL_SHAPES::Draw_Arrow(curve.particles.X(t),curve.particles.X(t)+velocity_scale*curve.particles.V(t));
        OpenGL_End();}

    glPopAttrib();

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SEGMENTED_CURVE_2D<T>::
Bounding_Box() const
{
    float xmin,xmax,ymin,ymax;
    xmin=xmax=curve.particles.X(1).x;
    ymin=ymax=curve.particles.X(1).y;
    for(int i=1; i<=curve.particles.array_collection->Size(); i++){
        xmin=min(xmin,(float)curve.particles.X(i).x);xmax=max(xmax,(float)curve.particles.X(i).x);
        ymin=min(ymin,(float)curve.particles.X(i).y);ymax=max(ymax,(float)curve.particles.X(i).y);}
    return World_Space_Box(RANGE<VECTOR<float,3> >(xmin,xmax,ymin,ymax,0,0));
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_SEGMENTED_CURVE_2D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION *selection = 0;
    if (buffer_size == 2)
    {
        if (buffer[0] == 1) selection = Get_Vertex_Selection(buffer[1]);
        else if (buffer[0] == 2) selection = Get_Segment_Selection(buffer[1]);
    }
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection; current_selection = 0;
    // Make a copy of selection
    if (selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D)
        current_selection=new OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T>(this,((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)selection)->index);
    else if (selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D)
        current_selection=new OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T>(this,((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)selection)->index);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const
{
    Print_Selection_Info(output_stream,selection,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection,MATRIX<T,3>* transform) const
{
    if (selection->type == OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D) {
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(index))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,curve.particles,index);}
    else if (selection->type == OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D) {
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)selection)->index;
        int node1,node2;curve.mesh.elements(index).Get(node1,node2);
        output_stream<<"Segment "<<index<<" ("<<node1<<", "<<node2<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(node1))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,curve.particles,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(node2))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,curve.particles,node2);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Clear_Highlight()
{
    delete current_selection; current_selection = 0;
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_SEGMENTED_CURVE_2D<T>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T>(this, index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_SEGMENTED_CURVE_2D<T>::
Get_Segment_Selection(int index)
{
    return new OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T>(this, index);
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T>
void OPENGL_SEGMENTED_CURVE_2D<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION::Draw_Vertices_For_Selection(curve.mesh,curve.particles);
}
//#####################################################################
// Function Draw_Segments_For_Selection
//#####################################################################
template<class T>
void OPENGL_SEGMENTED_CURVE_2D<T>::
Draw_Segments_For_Selection() const
{
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(OPENGL_PREFERENCES::selection_line_width);
    glPushName(0);
    for(int i=1;i<=curve.mesh.elements.m;i++){
        int node1,node2;curve.mesh.elements(i).Get(node1,node2);
        glLoadName(i);
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(curve.particles.X(node1),curve.particles.X(node2));
        OpenGL_End();
    }
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T>
RANGE<VECTOR<float,3> > OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const SEGMENTED_CURVE_2D<T> &curve=((OPENGL_SEGMENTED_CURVE_2D<T> *)object)->curve;
    RANGE<VECTOR<float,2> > box(VECTOR<float,2>(curve.particles.X(index)));
    return object->World_Space_Box(box);
}

template<class T>
RANGE<VECTOR<float,3> > OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const SEGMENTED_CURVE_2D<T> &curve=((OPENGL_SEGMENTED_CURVE_2D<T> *)object)->curve;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >::Bounding_Box(curve.particles.X.Subset(curve.mesh.elements(index))));
}
//#####################################################################
template class OPENGL_SEGMENTED_CURVE_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SEGMENTED_CURVE_2D<double>;
#endif
