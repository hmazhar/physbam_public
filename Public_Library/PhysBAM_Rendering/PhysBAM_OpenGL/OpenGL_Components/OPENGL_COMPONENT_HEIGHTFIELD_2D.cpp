//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_HEIGHTFIELD_2D.h>
using namespace PhysBAM;
//#####################################################################

template<class T,class RW> OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
OPENGL_COMPONENT_HEIGHTFIELD_2D(const GRID<TV> &grid_input, 
                                const std::string& height_filename_input,
                                const std::string& xz_filename_input,
                                const std::string& uv_filename_input,
                                int m_start_input,int m_end_input,int n_start_input,int n_end_input)
    :OPENGL_COMPONENT("Heightfield 2D"), 
    triangulated_surface(*TRIANGULATED_SURFACE<T>::Create()),
    opengl_triangulated_surface(triangulated_surface, false),
    vertical_offset(0), allow_smooth_shading(true), subdivide_surface(false),
    initial_grid(grid_input), grid(grid_input), height(grid.Domain_Indices()), 
    opengl_vector_field(vector_field, vector_locations, OPENGL_COLOR::Green(), 0.025, true, false, true),
    grid_filename(""), scale(1), displacement_scale(1), valid(false), draw_velocities(true), use_triangle_strip(false)
{
    if(m_start_input == 0 && m_end_input == 0) {
        domain=grid.Domain_Indices();
    } else {
        domain.min_corner.x = m_start_input; domain.max_corner.x = m_end_input; domain.min_corner.y = n_start_input; domain.max_corner.y = n_end_input;
    }
    counts=domain.Edge_Lengths()+1;
    PHYSBAM_ASSERT(counts.Min()>0);

    opengl_triangulated_surface.Set_Two_Sided();
    opengl_triangulated_surface.Set_Front_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0,0.6,(float)0.9)));
    opengl_triangulated_surface.Set_Back_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0.8,(float)0,0)));

    triangulated_surface.mesh.Initialize_Square_Mesh(counts.x,counts.y);
    for (int i=domain.min_corner.x;i<=domain.max_corner.x;i++) for(int j=domain.min_corner.y;j<=domain.max_corner.y;j++)
    {
        triangulated_surface.particles.array_collection->Add_Element();
    }

    height_filename=height_filename_input;
    if(xz_filename_input.length()){xz=new ARRAY<VECTOR<T,2> ,VECTOR<int,2> >;xz_filename=xz_filename_input;}
    else{xz=0;xz_filename="";}
    if(uv_filename_input.length()){
        uv = new ARRAY<VECTOR<T,2> ,VECTOR<int,2> >;
        vector_field.Resize(counts.Product());
        vector_locations.Resize(counts.Product());
        uv_filename=uv_filename_input;}
    else{uv=0;uv_filename="";}

    is_animation = FILE_UTILITIES::Is_Animated(height_filename);
    frame_loaded = -1;

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
~OPENGL_COMPONENT_HEIGHTFIELD_2D()
{
    delete &triangulated_surface.mesh;
    delete &triangulated_surface.particles;
    delete &triangulated_surface;
    delete xz;
}

template<class T,class RW> bool OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(is_animation?STRING_UTILITIES::string_sprintf(height_filename.c_str(),frame_input):height_filename);
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Display(const int in_color) const
{
    if(valid && draw)
    {
        if(use_triangle_strip)
        {
            opengl_triangulated_surface.front_material.Send_To_GL_Pipeline();
            VECTOR<T,3> v1,v2,v3,v4;
            for (int i = domain.min_corner.x; i <= domain.max_corner.x-1; i++)
            {
                OpenGL_Begin(GL_TRIANGLE_STRIP);
                v1=VECTOR<T,3>(grid.Axis_X(i+1,1), scale*(height(i+1,1)+vertical_offset), grid.Axis_X(1,2));
                v2=VECTOR<T,3>(grid.Axis_X(i,1), scale*(height(i,1)+vertical_offset), grid.Axis_X(1,2));
                OpenGL_Vertex(v1);
                OpenGL_Vertex(v2);
                for (int j = domain.min_corner.y+1; j <= domain.max_corner.y; j++)
                {
                    v3=VECTOR<T,3>(grid.Axis_X(i+1,1), scale*(height(i+1,j)+vertical_offset), grid.Axis_X(j,2));
                    v4=VECTOR<T,3>(grid.Axis_X(i,1), scale*(height(i,j)+vertical_offset), grid.Axis_X(j,2));

                    OpenGL_Normal(PLANE<T>(v1,v2,v3).Normal());
                    OpenGL_Vertex(v3);
                    OpenGL_Normal(PLANE<T>(v3,v2,v4).Normal());
                    OpenGL_Vertex(v4);
                    v1=v3;
                    v2=v4;
                }
                OpenGL_End();
            }
        }
        else
            opengl_triangulated_surface.Display(in_color);

        if(draw_velocities) opengl_vector_field.Display(in_color);
    }
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_triangulated_surface.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Turn_Smooth_Shading_On()
{
    if(allow_smooth_shading) opengl_triangulated_surface.Turn_Smooth_Shading_On();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Turn_Smooth_Shading_Off()
{
    opengl_triangulated_surface.Turn_Smooth_Shading_Off();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Reinitialize(bool force)
{
    if(draw){
        if(force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            bool success = true;
            valid = false;

            if(success && !grid_filename.empty()){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(grid_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,grid);
                    PHYSBAM_ASSERT(grid.counts.x==initial_grid.counts.x && grid.counts.y==initial_grid.counts.y);}
                else success=false;}

            if(success){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(height_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,height);
                    if(height.counts.x < counts.x || height.counts.y < counts.y) success = false;}
                else success=false;}

            if(success && xz){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(xz_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,*xz);
                    if(xz->counts.x < counts.x || xz->counts.y < counts.y) success = false;}
                else success=false;}

            if(success && draw_velocities && uv_filename.length()){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(uv_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    FILE_UTILITIES::Read_From_File<RW>(filename,*uv);
                    if(uv->counts.x < counts.x || uv->counts.y < counts.y)
                        success = false;
                    else{
                        int idx = 1;
                        for(int i=domain.min_corner.x;i<=domain.max_corner.x;i++)for(int j=domain.min_corner.y;j<=domain.max_corner.y;j++){
                            vector_field(idx) = VECTOR<T,3>((*uv)(i,j).x,0,(*uv)(i,j).y);
                            vector_locations(idx) = VECTOR<T,3>(grid.Axis_X(i,1), scale*(height(i,j)+vertical_offset), grid.Axis_X(j,2));
                            idx++;}}}
                else success=false;}

            if(success){
                Update_Surface();
                frame_loaded = frame;
                valid = true;}}}
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Update_Surface()
{
    if(use_triangle_strip) return;

    if(xz)
    {
        for (int i = domain.min_corner.x; i <= domain.max_corner.x; i++) for (int j = domain.min_corner.y; j <= domain.max_corner.y; j++)
        {
            triangulated_surface.particles.X(To_Linear_Index(i,j)) = 
                VECTOR<T,3>(grid.Axis_X(i,1) + displacement_scale*((*xz)(i,j).x-grid.Axis_X(i,1)), 
                             scale*(height(i,j)+vertical_offset), 
                             grid.Axis_X(j,2) + displacement_scale*((*xz)(i,j).y-grid.Axis_X(j,2)));
        }
    }
    else
    {
        if(subdivide_surface)
        {
            triangulated_surface.mesh.Initialize_Square_Mesh(grid.counts.x,grid.counts.y);
            triangulated_surface.particles.array_collection->Clean_Memory();
            triangulated_surface.particles.array_collection->Preallocate(grid.counts.x*grid.counts.y);
            for (int i = domain.min_corner.x; i <= domain.max_corner.x; i++) for (int j = domain.min_corner.y; j <= domain.max_corner.y; j++)
                triangulated_surface.particles.array_collection->Add_Element();
        }

        for (int i = domain.min_corner.x; i <= domain.max_corner.x; i++) for (int j = domain.min_corner.y; j <= domain.max_corner.y; j++)
        {
            triangulated_surface.particles.X(To_Linear_Index(i,j)) = 
                VECTOR<T,3>(grid.Axis_X(i,1), scale*(height(i,j)+vertical_offset), grid.Axis_X(j,2));
        }

        if(subdivide_surface) triangulated_surface.Loop_Subdivide();
    }
    if(opengl_triangulated_surface.Is_Smooth_Normals())
        opengl_triangulated_surface.Initialize_Vertex_Normals();
    else
        opengl_triangulated_surface.Delete_Vertex_Normals();
}

template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    if(use_triangle_strip) return 0;
    OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_2D<T> *new_selection = 0;
    OPENGL_SELECTION *selection = opengl_triangulated_surface.Get_Selection(buffer,buffer_size);
    if(selection)
    {
        if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX)
        {
            int index = ((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T> *)selection)->index;
            new_selection = new OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_2D<T>(this);
            new_selection->index = From_Linear_Index(index);
        }
        delete selection;
    }
    return new_selection;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type != OPENGL_SELECTION::COMPONENT_HEIGHTFIELD_2D) return;
    OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_2D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_2D<T>*)selection;
    int vertex_index = To_Linear_Index(real_selection->index.x, real_selection->index.y);
    OPENGL_SELECTION *surface_selection = opengl_triangulated_surface.Get_Vertex_Selection(vertex_index);
    opengl_triangulated_surface.Highlight_Selection(surface_selection);
    delete surface_selection; // Highlight_Selection made a copy of it
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Clear_Highlight()
{
    opengl_triangulated_surface.Clear_Highlight();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Set_Scale(T scale_input)
{
    scale = scale_input;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Use_Triangle_Strip(bool use_triangle_strip_input)
{
    use_triangle_strip=use_triangle_strip_input;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Increase_Scale()
{
    scale *= 1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Decrease_Scale()
{
    scale *= 1/1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Increase_Displacement_Scale()
{
    displacement_scale *= 1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Decrease_Displacement_Scale()
{
    displacement_scale *= 1/1.1;
    if(valid) Update_Surface();
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Increase_Velocity_Scale()
{
    opengl_vector_field.size *= 1.1;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Decrease_Velocity_Scale()
{
    opengl_vector_field.size *= 1/1.1;
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Toggle_Draw_Velocities()
{
    draw_velocities = !draw_velocities;
    Reinitialize(true);
}

template<class T,class RW> void OPENGL_COMPONENT_HEIGHTFIELD_2D<T,RW>::
Toggle_Subdivision()
{
    subdivide_surface = !subdivide_surface;
    Reinitialize(true);
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_HEIGHTFIELD_2D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<float,3> >::Centered_Box();
}

template class OPENGL_COMPONENT_HEIGHTFIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_HEIGHTFIELD_2D<double,double>;
#endif
