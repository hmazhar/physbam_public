//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_HEIGHTFIELD_1D.h>
using namespace PhysBAM;
//#####################################################################

template<class T,class RW> OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
OPENGL_COMPONENT_HEIGHTFIELD_1D(const GRID<TV> &grid_input, 
                                const std::string& height_filename_input,
                                const std::string& x_filename_input,
                                const std::string& ground_filename_input,
                                const std::string& u_filename_input)
    : OPENGL_COMPONENT("Heightfield 1D"), grid(grid_input), opengl_vector_field(vector_field,vector_locations), scale(1), displacement_scale(1), valid(false), 
      draw_velocities(true), draw_points(true), selected_index(0)
{
    height_filename=height_filename_input;
    if(x_filename_input.length()){x=new ARRAY<T,TV_INT>;x_filename=x_filename_input;}
    else{x=0;x_filename="";}
    if(ground_filename_input.length()){ground=new ARRAY<T,TV_INT>;ground_filename=ground_filename_input;}
    else{ground=0;ground_filename="";}
    if(u_filename_input.length()){
        u = new ARRAY<T,TV_INT>;
        vector_field.Resize(grid.counts.x);
        vector_locations.Resize(grid.counts.x);
        u_filename=u_filename_input;}
    else{u=0;u_filename="";}
    is_animation=height_filename.find("%d")!=std::string::npos;
    frame_loaded = -1;

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
~OPENGL_COMPONENT_HEIGHTFIELD_1D()
{
    delete x;
    delete ground;
    delete u;
}

template<class T,class RW> bool OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(is_animation?STRING_UTILITIES::string_sprintf(height_filename.c_str(),frame_input):height_filename);
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Display(const int in_color) const
{
    OPENGL_COLOR water_color(0,0,1);
    OPENGL_COLOR ground_color = OPENGL_COLOR::Gray(0.8);
    OPENGL_COLOR displaced_water(1,0,0);
    OPENGL_COLOR points_color(0,1,1);
    OPENGL_COLOR selected_point_color = OPENGL_COLOR::Yellow();
    GLfloat point_size = 3.0;

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    if (valid && draw)
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);

        if(mode == GL_SELECT)
        {
            glPushName(0);
            glPushAttrib(GL_POINT_BIT);
            glPointSize(point_size);
            points_color.Send_To_GL_Pipeline();
            for (int i = 1; i <= grid.counts.x; i++)
            {
                glLoadName(i);
                OpenGL_Begin(GL_POINTS);
                OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(i,1), scale*height(i)));
                OpenGL_End();
            }
            glPopAttrib();
            glPopName();
        }
        else
        {
            if (x)
            {
                displaced_water.Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINE_STRIP);
                if (displacement_scale == 1)
                {
                    for (int i = 1; i <= x->counts.x; i++)
                        OpenGL_Vertex(VECTOR<T,2>((*x)(i), scale*height(i)));        
                }
                else
                {
                    for (int i = 1; i <= x->counts.x; i++)
                        OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(i,1) + displacement_scale*((*x)(i)-grid.Axis_X(i,1)), 
                                                   scale*height(i)));        
                }
                OpenGL_End();
            }

            // Water
            water_color.Send_To_GL_Pipeline();
            OpenGL_Begin(GL_LINE_STRIP);
            for (int i = 1; i <= grid.counts.x; i++)
                OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(i,1), scale*height(i)));
            OpenGL_End();

            if (draw_points)
            {
                glPushAttrib(GL_POINT_BIT);
                glPointSize(point_size);
                points_color.Send_To_GL_Pipeline();
                OpenGL_Begin(GL_POINTS);
                for (int i = 1; i <= grid.counts.x; i++)
                    OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(i,1), scale*height(i)));
                OpenGL_End();

                for(int i=1;i<=grid.counts.x;i++)
                    OpenGL_String(VECTOR<T,2>(grid.Axis_X(i,1),scale*height(i)),STRING_UTILITIES::string_sprintf("%d",i));

                if(selected_index){
                    selected_point_color.Send_To_GL_Pipeline();
                    int i=selected_index;
                    OpenGL_Begin(GL_POINTS);
                    OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(i,1), scale*height(i)));
                    OpenGL_End();
                    OpenGL_String(VECTOR<T,2>(grid.Axis_X(i,1),scale*height(i)),STRING_UTILITIES::string_sprintf("%d",i));}

                glPopAttrib();
            }

            // Ground
            ground_color.Send_To_GL_Pipeline();
            if (ground)
            {
                OpenGL_Begin(GL_LINE_STRIP);
                for (int i = 1; i <= grid.counts.x; i++)
                    OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(i,1), scale*(*ground)(i)));        
                OpenGL_End();
            }
            else
            {
                OpenGL_Begin(GL_LINE_STRIP);
                OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(1,1), 0));
                OpenGL_Vertex(VECTOR<T,2>(grid.Axis_X(grid.counts.x,1), 0));
                OpenGL_End();
            }
        }

        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);

        if (mode != GL_SELECT && draw_velocities) opengl_vector_field.Display(in_color);
    }
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Bounding_Box() const
{
    if (valid && draw)
    {
        T min_height=min((T)0,height.Min());
        T max_height=max((T)0,height.Max());

        T min_x = grid.domain.min_corner.x, max_x = grid.domain.max_corner.x;
        if (x) { 
            min_x=min(min_x,x->Min()); 
            max_x=max(max_x,x->Max());
        }

        return RANGE<VECTOR<float,3> >(min_x, max_x, min_height, max_height, 0, 0);
    }
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Reinitialize(bool force)
{
    if(draw){
        if(force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            bool success = true;
            valid = false;

            if (success){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(height_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,height);
                    if(height.counts.x != grid.counts.x) success=false;}
                else success=false;}

            if(success && x){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(x_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,*x);
                    if(height.counts.x != x->counts.x) success=false;}
                else success=false;}
    
            if(success && ground){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(ground_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,*ground);
                    if (height.counts.x != ground->counts.x) success = false;
                    else height += (*ground);}
                else success=false;}

            if(success && draw_velocities && u_filename.length()){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(u_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,*u);
                    if (height.counts.x != u->counts.x) success = false;
                    else for(int i=1;i<=grid.counts.x;i++){
                            vector_field(i) = VECTOR<T,2>((*u)(i),0);
                            vector_locations(i) = VECTOR<T,2>(grid.Axis_X(i,1), scale*height(i));}}
                else success=false;}

            if(success){
                frame_loaded = frame;
                valid = true;}}}
}

template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    if (buffer_size == 1)
    {
        OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T> *selection = new OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T>(this);
        selection->index = buffer[0];
        return selection;
    }
    else return 0;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if (selection->type != OPENGL_SELECTION::COMPONENT_HEIGHTFIELD_1D) return;
    OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T>*)selection;
    selected_index = real_selection->index;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Clear_Highlight()
{
    selected_index = 0;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Set_Scale(T scale_input)
{
    scale = scale_input;
    Reinitialize(true); // To recompute velocity vector positions correctly
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Increase_Scale()
{
    scale *= 1.1;
    Reinitialize(true); // To recompute velocity vector positions correctly
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Decrease_Scale()
{
    scale *= 1/1.1;
    Reinitialize(true); // To recompute velcity vector positions correctly
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Increase_Displacement_Scale()
{
    displacement_scale *= 1.1;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Decrease_Displacement_Scale()
{
    displacement_scale *= 1/1.1;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Increase_Velocity_Scale()
{
    opengl_vector_field.size *= 1.1;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Decrease_Velocity_Scale()
{
    opengl_vector_field.size *= 1/1.1;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Toggle_Draw_Velocities()
{
    draw_velocities = !draw_velocities;
    Reinitialize(true);
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_1D<T,RW>::
Toggle_Draw_Points()
{
    draw_points = !draw_points;
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_1D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Empty_Box();
}

template class OPENGL_COMPONENT_HEIGHTFIELD_1D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_HEIGHTFIELD_1D<double,double>;
#endif
