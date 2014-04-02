//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Geoffrey Irving, Michael Lentine, Neil Molino, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
using namespace PhysBAM;
OPENGL_OBJECT::
OPENGL_OBJECT()
    :frame(&default_frame),selectable(false),visible(true),show_name(true),slice(0)
{
}
OPENGL_OBJECT::
~OPENGL_OBJECT()
{
}
void OPENGL_OBJECT::
Display(const int in_color) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
bool OPENGL_OBJECT::
Use_Bounding_Box() const
{
    return true;
}
RANGE<VECTOR<float,3> > OPENGL_OBJECT::
Bounding_Box() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
bool OPENGL_OBJECT::
Is_Transparent() const
{
    return false;  // this is usually what we want
}
void OPENGL_OBJECT::
Turn_Smooth_Shading_On()
{
}
void OPENGL_OBJECT::
Turn_Smooth_Shading_Off()
{
}
OPENGL_SELECTION *OPENGL_OBJECT::
Get_Selection(GLuint *buffer,int buffer_size)
{
    return 0;
}
void OPENGL_OBJECT::
Set_Selection(OPENGL_SELECTION *selection)
{
}
void OPENGL_OBJECT::
Highlight_Selection(OPENGL_SELECTION *selection)
{
}
void OPENGL_OBJECT::
Clear_Highlight()
{
}
void OPENGL_OBJECT::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
}
RANGE<VECTOR<float,3> > OPENGL_OBJECT::
Selection_Bounding_Box(OPENGL_SELECTION *selection) const
{
    return RANGE<VECTOR<float,3> >::Centered_Box();
}
OPENGL_SELECTION* OPENGL_OBJECT::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    return 0;
}
void OPENGL_OBJECT::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;Slice_Has_Changed();
}
void OPENGL_OBJECT::
Slice_Has_Changed()
{
}
