//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D(const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT("Face Scalar Field 1D"),opengl_face_scalar_field(grid_input,*new ARRAY<T2,FACE_INDEX<1> >,point_color,line_color),
      values_filename(values_filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(values_filename);
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
~OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D()
{
    delete &opengl_face_scalar_field.face_values;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2,class RW> bool OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(values_filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_face_scalar_field.Display(in_color);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_face_scalar_field.Print_Selection_Info(stream,selection);}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_face_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Scale
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Scale(const T scale)
{
    opengl_face_scalar_field.Scale(scale);
}
//#####################################################################
// Function Increase_Scale
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Increase_Scale()
{
    Scale((T)2);
}
//#####################################################################
// Function Decrease_Scale 
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Decrease_Scale()
{
    Scale((T).5);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2,RW>::
Reinitialize()
{
    if(draw && ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))){
        valid=false;

        std::string filename=FILE_UTILITIES::Get_Frame_Filename(values_filename,frame);
        if(FILE_UTILITIES::File_Exists(filename))
            FILE_UTILITIES::Read_From_File<RW>(filename,opengl_face_scalar_field.face_values);
        else return;
        frame_loaded = frame;
        valid = true;}
}
//#####################################################################
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<float,int,float>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<float,bool,float>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<double,int,double>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<double,bool,double>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<double,double,double>;
#endif
