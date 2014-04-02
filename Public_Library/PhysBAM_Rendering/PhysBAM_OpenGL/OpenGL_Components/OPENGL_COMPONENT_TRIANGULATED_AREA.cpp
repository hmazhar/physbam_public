//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
OPENGL_COMPONENT_TRIANGULATED_AREA(const std::string &filename)
    : OPENGL_COMPONENT("Triangulated Surface"),color_map(0), 
      triangulated_area(*TRIANGULATED_AREA<T>::Create()),
      opengl_triangulated_area(triangulated_area),
      filename(filename),color_map_filename(0),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
OPENGL_COMPONENT_TRIANGULATED_AREA(const std::string &filename,const std::string &color_map_filename_input)
    : OPENGL_COMPONENT("Triangulated Surface"),color_map(new ARRAY<OPENGL_COLOR >), 
      triangulated_area(*TRIANGULATED_AREA<T>::Create()),
      opengl_triangulated_area(triangulated_area,false,OPENGL_COLOR::Red(),OPENGL_COLOR::Black()),
      filename(filename),color_map_filename(&color_map_filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    opengl_triangulated_area.Set_Color_Map(color_map);
    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
~OPENGL_COMPONENT_TRIANGULATED_AREA()
{
    delete &triangulated_area.mesh;
    delete &triangulated_area.particles;
    delete &triangulated_area;
}

template<class T,class RW> bool OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_triangulated_area.Display(in_color);
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_triangulated_area.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA<T,RW>::
Reinitialize()
{
    if (draw)
    {
        if ((is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(filename, frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,triangulated_area);
            else
                return;
            if(color_map) {
                std::string tmp_color_map_filename = FILE_UTILITIES::Get_Frame_Filename(*color_map_filename, frame);
                //if (FILE_UTILITIES::File_Exists(tmp_filename))
                if (FILE_UTILITIES::File_Exists(tmp_color_map_filename))
                    FILE_UTILITIES::Read_From_File<RW>(tmp_color_map_filename,*color_map);
                else
                    return;
            }            
            frame_loaded = frame;
            valid = true;
        }
    }
}
template class OPENGL_COMPONENT_TRIANGULATED_AREA<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_TRIANGULATED_AREA<double,double>;
#endif
