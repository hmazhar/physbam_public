//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DIAGNOSTICS.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// OPENGL_COMPONENT_DIAGNOSTICS
//#####################################################################
OPENGL_COMPONENT_DIAGNOSTICS::
OPENGL_COMPONENT_DIAGNOSTICS(const std::string& filename_input)
    :filename(filename_input),frame_loaded(INT_MIN),valid(false)
{
}
//#####################################################################
// ~OPENGL_COMPONENT_DIAGNOSTICS
//#####################################################################
OPENGL_COMPONENT_DIAGNOSTICS::
~OPENGL_COMPONENT_DIAGNOSTICS()
{
}
//#####################################################################
// Valid_Frame
//#####################################################################
bool OPENGL_COMPONENT_DIAGNOSTICS::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
void OPENGL_COMPONENT_DIAGNOSTICS::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
void OPENGL_COMPONENT_DIAGNOSTICS::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    output_stream<<component_name<<":"<<std::endl;
    for(int i=1;i<=lines.m;i++) output_stream<<"   "<<lines(i)<<std::endl;
}
//#####################################################################
// Reinitialize
//#####################################################################
void OPENGL_COMPONENT_DIAGNOSTICS::
Reinitialize()
{
    if(draw){
      if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded==INT_MIN)){
            lines.Remove_All();
            valid=false;
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(filename, frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename)){
                std::istream* input=FILE_UTILITIES::Safe_Open_Input(tmp_filename,false);
                std::string line;while(std::getline(*input,line)) lines.Append(line);
                delete input;frame_loaded=frame;valid=true;}}}
}
//#####################################################################
// Use_Bounding_Box
//#####################################################################
bool OPENGL_COMPONENT_DIAGNOSTICS::
Use_Bounding_Box() const
{
    return false;
}
//#####################################################################
// Display
//#####################################################################
void OPENGL_COMPONENT_DIAGNOSTICS::
Display(const int in_color) const
{
}
//#####################################################################
