//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Reinitialize
//##################################################################### 
template<class T,class RW> void OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T,RW>::
Reinitialize(bool force_load_even_if_not_drawn)
{
    if(!draw&&!force_load_even_if_not_drawn)return;
    if(!valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
        valid=false;std::string tmp_filename=FILE_UTILITIES::Get_Frame_Filename(field_filename,frame);
        if(FILE_UTILITIES::File_Exists(tmp_filename))FILE_UTILITIES::Read_From_File<RW>(tmp_filename,field);else return;
        opengl_symmetric_matrix_field.Update();
        frame_loaded=frame;valid=true;}
}
//##################################################################### 
template class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<double,double>;
#endif
