//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_1D
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
OPENGL_COMPONENT_LEVELSET_1D(GRID<TV> &grid,const std::string& levelset_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT("Levelset 1D"),levelset_filename(levelset_filename_input),opengl_levelset(0)
{
    is_animation=FILE_UTILITIES::Is_Animated(levelset_filename);
    opengl_levelset=new OPENGL_LEVELSET_1D<T>(*(new LEVELSET_1D<T>(grid,*(new T_ARRAYS_SCALAR))),point_color,line_color);
    Reinitialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
~OPENGL_COMPONENT_LEVELSET_1D()
{
    delete &opengl_levelset->levelset;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(levelset_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Display(const int in_color) const
{
    opengl_levelset->Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_levelset->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Reinitialize()
{
    if(draw && ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))){
        valid=false;
        std::string filename=FILE_UTILITIES::Get_Frame_Filename(levelset_filename,frame);
        if(FILE_UTILITIES::File_Exists(filename)){
            FILE_UTILITIES::Read_From_File<RW>(filename,opengl_levelset->levelset);
            frame_loaded=frame;valid=true;}}
}
//##################################################################### 
template class OPENGL_COMPONENT_LEVELSET_1D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_LEVELSET_1D<float,double>;
template class OPENGL_COMPONENT_LEVELSET_1D<double,double>;
#endif
