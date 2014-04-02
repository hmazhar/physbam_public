//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################

template<class T,class T2,class RW> OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D(const GRID<TV> &grid_input, const std::string &values_filename_input, OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT("Face Scalar Field 2D"), opengl_scalar_field(grid_input,*new ARRAY<T2,FACE_INDEX<2> >,color_map_input),
    values_filename(values_filename_input), x_face_values_filename(""), y_face_values_filename(""),  frame_loaded(-1), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(values_filename);
    Reinitialize();
}

template<class T,class T2,class RW> OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D(const GRID<TV> &grid_input, const std::string &x_face_values_filename_input, const std::string &y_face_values_filename_input,
                                      OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT("Face Scalar Field 2D"), opengl_scalar_field(grid_input,*new ARRAY<T2,FACE_INDEX<2> >,color_map_input),
    values_filename(), x_face_values_filename(x_face_values_filename_input), y_face_values_filename(y_face_values_filename_input),
    frame_loaded(-1), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(x_face_values_filename);
    Reinitialize();
}

template<class T,class T2,class RW> OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
~OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D()
{
    delete &opengl_scalar_field.x_face_values;
    delete &opengl_scalar_field.y_face_values;
}

template<class T,class T2,class RW> bool OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    if (!values_filename.empty())
        return FILE_UTILITIES::Frame_File_Exists(values_filename, frame_input);
    else
        return FILE_UTILITIES::Frame_File_Exists(x_face_values_filename, frame_input) && 
               FILE_UTILITIES::Frame_File_Exists(y_face_values_filename, frame_input);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_scalar_field.Display(in_color);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_scalar_field.Print_Selection_Info(stream,selection);}
}

template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2,RW>::
Reinitialize()
{
    if (draw)
    {
        if ((is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string filename;
            if(!values_filename.empty()){
                filename=FILE_UTILITIES::Get_Frame_Filename(values_filename, frame);
                if(FILE_UTILITIES::File_Exists(filename))
                    FILE_UTILITIES::Read_From_File<RW>(filename,opengl_scalar_field.face_values);
                else return;}
            else{
                filename=FILE_UTILITIES::Get_Frame_Filename(x_face_values_filename, frame);
                if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename,opengl_scalar_field.x_face_values);
                else return;
                filename=FILE_UTILITIES::Get_Frame_Filename(y_face_values_filename, frame);
                if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename,opengl_scalar_field.y_face_values);
                else return;
            }

            opengl_scalar_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<float,int,float>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<float,bool,float>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<double,int,double>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<double,bool,double>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<double,double,double>;
#endif
