//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
OPENGL_COMPONENT_TRIANGULATED_SURFACE(const std::string &filename, bool use_display_list)
    : OPENGL_COMPONENT("Triangulated Surface"), 
      triangulated_surface(*TRIANGULATED_SURFACE<T>::Create()),
      opengl_triangulated_surface(triangulated_surface, false,
                                  OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Red()),
                                  OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Blue())),
      filename(filename), frame_loaded(-1), valid(false), use_display_list(use_display_list)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
~OPENGL_COMPONENT_TRIANGULATED_SURFACE()
{
    delete &triangulated_surface.mesh;
    delete &triangulated_surface.particles;
    delete &triangulated_surface;
}

template<class T,class RW> bool OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
Display(const int in_color) const
{
    if (valid && draw) 
    {
        if (slice && slice->Is_Slice_Mode()) {
            glPushAttrib(GL_ENABLE_BIT);
            slice->Enable_Clip_Planes();
        }
        opengl_triangulated_surface.Display(in_color);
        if (slice && slice->Is_Slice_Mode()) {
            glPopAttrib();
        }
    }
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_triangulated_surface.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,RW>::
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
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,triangulated_surface);
            else
                return;

            opengl_triangulated_surface.name=tmp_filename;
            frame_loaded = frame;
            valid = true;
        }
    }
}
template class OPENGL_COMPONENT_TRIANGULATED_SURFACE<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_TRIANGULATED_SURFACE<double,double>;
#endif
