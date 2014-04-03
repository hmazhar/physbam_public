//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD(TRIANGULATED_AREA<T>& triangulated_area,const std::string &vector_field_filename)
    :OPENGL_COMPONENT("Triangulated Area Based Vector Field 2D"),opengl_vector_field(triangulated_area,*new ARRAY<VECTOR<T,2> >), 
    vector_field_filename(vector_field_filename),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(vector_field_filename);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
~OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD()
{
    delete &opengl_vector_field.V;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(vector_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_vector_field.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_vector_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Reinitialize(bool force_load_even_if_not_drawn)
{
    if (draw||force_load_even_if_not_drawn)
    {
        if (!valid ||
            (is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(vector_field_filename, frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_vector_field.V);
            else
                return;

            opengl_vector_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Increase_Vector_Size()
{
    opengl_vector_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Decrease_Vector_Size()
{
    opengl_vector_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>::
Toggle_Arrowhead()
{
    opengl_vector_field.draw_arrowhead = !opengl_vector_field.draw_arrowhead;
}
template class OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<double,double>;
#endif
