//#####################################################################
// Copyright 2006-2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SCALAR_FIELD_1D
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_1D.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
OPENGL_COMPONENT_SCALAR_FIELD_1D(const GRID<TV> &grid,const std::string &scalar_field_filename,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT("Scalar Field 1D"),scalar_field_filename(scalar_field_filename),frame_loaded(INT_MIN),valid(false),opengl_scalar_field(grid,*new ARRAY<T2,VECTOR<int,1> >,point_color,line_color)
{
    is_animation=FILE_UTILITIES::Is_Animated(scalar_field_filename);
}
//#####################################################################
// Function ~OPENGL_COMPONENT_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
~OPENGL_COMPONENT_SCALAR_FIELD_1D()
{
    delete &opengl_scalar_field.values;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2,class RW> bool OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Display(const int in_color) const
{
    if(valid && draw) opengl_scalar_field.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Reinitialize()
{
    if(draw && ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))){
        valid=false;
        std::string filename=FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename,frame);
        if(FILE_UTILITIES::File_Exists(filename)){
            FILE_UTILITIES::Read_From_File<RW>(filename,opengl_scalar_field.values);
            frame_loaded=frame;valid=true;}}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Scale(const T scale)
{
    opengl_scalar_field.Scale(scale);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Increase_Scale()
{
    Scale((T)2);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2,RW>::
Decrease_Scale()
{
    Scale((T).5);
}
//##################################################################### 
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<float,float,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<float,float,double>;
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<float,bool,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<float,bool,double>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<double,bool,double>;
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<double,double,double>;
#endif
