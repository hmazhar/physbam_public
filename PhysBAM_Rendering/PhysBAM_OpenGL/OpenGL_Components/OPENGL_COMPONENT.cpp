//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
using namespace PhysBAM;
OPENGL_COMPONENT::
OPENGL_COMPONENT(const std::string &name)
 : frame(0), draw(true), is_animation(false), component_name(name)
{
}
OPENGL_COMPONENT::
~OPENGL_COMPONENT()
{
}
bool OPENGL_COMPONENT::
Valid_Frame(int frame_input) const
{
    return false;
}
bool OPENGL_COMPONENT::
Is_Up_To_Date(int frame) const
{
    return true;
}
void OPENGL_COMPONENT::
Set_Frame(int frame_input)
{
    frame = frame_input;
}
void OPENGL_COMPONENT::
Set_Draw(bool draw_input)
{
    draw = draw_input;
}
void OPENGL_COMPONENT::
Draw_All_Objects()
{
}
bool OPENGL_COMPONENT::
Use_Bounding_Box() const
{
    return draw;
}
void OPENGL_COMPONENT::
Next_Frame()
{
    Set_Frame(frame+1);
}
void OPENGL_COMPONENT::
Prev_Frame()
{
    Set_Frame(frame-1);
}
