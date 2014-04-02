//#####################################################################
// Copyright 2004, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TRIANGULATED_AREA<T>::
OPENGL_TRIANGULATED_AREA(TRIANGULATED_AREA<T>& triangulated_area_input,const bool draw_vertices_input,const OPENGL_COLOR& vertex_color_input,
    const OPENGL_COLOR& segment_color_input,const OPENGL_COLOR& triangle_color_input,ARRAY<OPENGL_COLOR>* color_map_input)
    :triangulated_area(triangulated_area_input),vertex_color(vertex_color_input),segment_color(segment_color_input),triangle_color(triangle_color_input),
    velocity_color(OPENGL_COLOR::Yellow()),current_selection(0),color_map(color_map_input),draw_vertices(draw_vertices_input),draw_velocities(false),velocity_scale(0.025)
{}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    if(mode == GL_SELECT){
        glPushName(1);
        Draw_Vertices_For_Selection();
        if(triangulated_area.mesh.segment_mesh){
            glLoadName(2);
            Draw_Segments_For_Selection();}
        glLoadName(3);
        Draw_Triangles_For_Selection();
        glPopName();}
    else{
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_POLYGON_BIT);
        glDisable(GL_LIGHTING);glEnable(GL_POLYGON_OFFSET_FILL);glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(0,2);
        triangle_color.Send_To_GL_Pipeline();
        Draw_Triangles();
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        glPolygonOffset(0,1);
        segment_color.Send_To_GL_Pipeline();
        Draw_Triangles(false);
        Draw_Vertices();
        glPopAttrib();}

    if(mode != GL_SELECT){
        if(current_selection){
            glPushAttrib(GL_ENABLE_BIT);
            glDisable(GL_DEPTH_TEST);
            if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX){
                int index=((OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T> *)current_selection)->index;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(triangulated_area.particles.X(index));}
            else if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT && triangulated_area.mesh.segment_mesh){
                int index=((OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T> *)current_selection)->index;
                int node1,node2;triangulated_area.mesh.segment_mesh->elements(index).Get(node1,node2);
                OPENGL_SELECTION::Draw_Highlighted_Segment(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2));}
            else if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE){
                int index=((OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T> *)current_selection)->index;
                int node1,node2,node3;triangulated_area.mesh.elements(index).Get(node1,node2,node3);
                OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2),triangulated_area.particles.X(node3));}
            glPopAttrib();
        }
    }
    if (draw_velocities && triangulated_area.particles.store_velocity){
        glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int t=1;t<=triangulated_area.particles.array_collection->Size();t++)
            OPENGL_SHAPES::Draw_Arrow(triangulated_area.particles.X(t),triangulated_area.particles.X(t)+velocity_scale*triangulated_area.particles.V(t));
        OpenGL_End();
        glPopAttrib();}

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_TRIANGULATED_AREA<T>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box=RANGE<VECTOR<float,3> >::Empty_Box();
    for(int i=1;i<=triangulated_area.particles.array_collection->Size();i++)box.Enlarge_To_Include_Point(World_Space_Point(VECTOR<float,2>(triangulated_area.particles.X(i))));
    return box;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TRIANGULATED_AREA<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size == 2){
        if(buffer[0] == 1) selection = Get_Vertex_Selection(buffer[1]);
        else if(buffer[0] == 2) selection = Get_Segment_Selection(buffer[1]);
        else if(buffer[0] == 3) selection = Get_Triangle_Selection(buffer[1]);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Highlight_Selection(OPENGL_SELECTION* selection)
{
    delete current_selection; current_selection = 0;
    // Make a copy of selection
    if(selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX)
        current_selection=new OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T>(this,((OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T>*)selection)->index);
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT)
        current_selection=new OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T>(this,((OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T>*)selection)->index);
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE)
        current_selection=new OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T>(this,((OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T>*)selection)->index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection) const
{
    Print_Selection_Info(output_stream,selection,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection,MATRIX<T,3>* transform) const
{
    if(selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX){
        int index=((OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T>*)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(index))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,triangulated_area.particles,index);}
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT){
        PHYSBAM_ASSERT(triangulated_area.mesh.segment_mesh);
        int index=((OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T>*)selection)->index;
        int node1,node2;triangulated_area.mesh.segment_mesh->elements(index).Get(node1,node2);
        output_stream<<"Segment "<<index<<" ("<<node1<<", "<<node2<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node1))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,triangulated_area.particles,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node2))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,triangulated_area.particles,node2);}
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE){
        int index=((OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T>*)selection)->index;
        int node1,node2,node3;triangulated_area.mesh.elements(index).Get(node1,node2,node3);
        output_stream<<"Triangle "<<index<<" ("<<node1<<", "<<node2<<", "<<node3<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node1))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,triangulated_area.particles,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node2))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,triangulated_area.particles,node2);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node3<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node3))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,2> >,T>::Print(output_stream,triangulated_area.particles,node3);}
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TRIANGULATED_AREA<T>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T>(this,index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TRIANGULATED_AREA<T>::
Get_Segment_Selection(int index)
{
    return new OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T>(this,index);
}
//#####################################################################
// Function Get_Triangle_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TRIANGULATED_AREA<T>::
Get_Triangle_Selection(int index)
{
    return new OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T>(this,index);
}
//#####################################################################
// Function Draw_Vertices
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Vertices() const
{
    if(draw_vertices){
        vertex_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_POINTS);
        if(!triangulated_area.mesh.neighbor_nodes) triangulated_area.mesh.Initialize_Neighbor_Nodes();
        for(int i=1;i<=triangulated_area.particles.array_collection->Size();i++) if((*triangulated_area.mesh.neighbor_nodes)(i).m)
            OpenGL_Vertex(triangulated_area.particles.X(i));
        OpenGL_End();}
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION::Draw_Vertices_For_Selection(triangulated_area.mesh,triangulated_area.particles);
}
//#####################################################################
// Function Draw_Segments
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Segments() const
{
    PHYSBAM_ASSERT(triangulated_area.mesh.segment_mesh);
    segment_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINES);
    for(int i=1;i<=triangulated_area.mesh.segment_mesh->elements.m;i++){
        int node1,node2;triangulated_area.mesh.segment_mesh->elements(i).Get(node1,node2);
        OpenGL_Line(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2));}
    OpenGL_End();
}
//#####################################################################
// Function Draw_Segments_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Segments_For_Selection() const
{
    PHYSBAM_ASSERT(triangulated_area.mesh.segment_mesh);
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(OPENGL_PREFERENCES::selection_line_width);
    glPushName(0);
    for(int i=1;i<=triangulated_area.mesh.segment_mesh->elements.m;i++){
        int node1,node2;
        triangulated_area.mesh.segment_mesh->elements(i).Get(node1,node2);
        glLoadName(i);
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2));
        OpenGL_End();}
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Triangles
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Triangles(const bool use_color_map) const
{
    OpenGL_Begin(GL_TRIANGLES);
    for(int i=1;i<=triangulated_area.mesh.elements.m;i++){
        if(color_map && use_color_map) (*color_map)(i).Send_To_GL_Pipeline();
        int node1,node2,node3;triangulated_area.mesh.elements(i).Get(node1,node2,node3);
        OpenGL_Triangle(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2),triangulated_area.particles.X(node3));}
    OpenGL_End();
}
//#####################################################################
// Function Draw_Triangles_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Triangles_For_Selection() const
{
    glPushName(0);
    for(int i=1;i<=triangulated_area.mesh.elements.m;i++){
        int node1,node2,node3;triangulated_area.mesh.elements(i).Get(node1,node2,node3);
        glLoadName(i);
        OpenGL_Begin(GL_TRIANGLES);
        OpenGL_Triangle(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2),triangulated_area.particles.X(node3));
        OpenGL_End();}
    glPopName();
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const TRIANGULATED_AREA<T>& triangulated_area=((OPENGL_TRIANGULATED_AREA<T>*)object)->triangulated_area;
    RANGE<VECTOR<T,2> > box(triangulated_area.particles.X(index));
    return object->World_Space_Box(RANGE<VECTOR<float,2> >(box));
}
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const TRIANGULATED_AREA<T>& triangulated_area=((OPENGL_TRIANGULATED_AREA<T>*)object)->triangulated_area;
    PHYSBAM_ASSERT(triangulated_area.mesh.segment_mesh);
    return object->World_Space_Box(RANGE<VECTOR<T,2> >::Bounding_Box(triangulated_area.particles.X.Subset(triangulated_area.mesh.segment_mesh->elements(index))));
}
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const TRIANGULATED_AREA<T>& triangulated_area=((OPENGL_TRIANGULATED_AREA<T>*)object)->triangulated_area;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >::Bounding_Box(triangulated_area.particles.X.Subset(triangulated_area.mesh.elements(index))));
}
//#####################################################################
template class OPENGL_TRIANGULATED_AREA<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_TRIANGULATED_AREA<double>;
#endif
