//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
using namespace PhysBAM;

// Default (non-smooth rendering, thin lines)
bool OPENGL_PREFERENCES::smooth_points=false;
bool OPENGL_PREFERENCES::smooth_lines=false;
float OPENGL_PREFERENCES::line_width=1;
float OPENGL_PREFERENCES::point_size=2;
float OPENGL_PREFERENCES::selection_point_size=10;
float OPENGL_PREFERENCES::selection_line_width=3;
float OPENGL_PREFERENCES::highlighted_point_size=10;
float OPENGL_PREFERENCES::highlighted_line_width=3;

const float OPENGL_PREFERENCES::arrowhead_size=0.3;
const float OPENGL_PREFERENCES::arrowhead_angle=0.3; // radians
const OPENGL_COLOR OPENGL_PREFERENCES::selection_highlight_color=OPENGL_COLOR::Yellow();

// Preferences for smooth rendering
void OPENGL_PREFERENCES::Set_Smooth_Defaults()
{
    OPENGL_PREFERENCES::smooth_points=true;
    OPENGL_PREFERENCES::smooth_lines=true;
    OPENGL_PREFERENCES::line_width=2;
    OPENGL_PREFERENCES::point_size=3;
    OPENGL_PREFERENCES::selection_point_size=10;
    OPENGL_PREFERENCES::selection_line_width=3;
    OPENGL_PREFERENCES::highlighted_point_size=10;
    OPENGL_PREFERENCES::highlighted_line_width=6;
}
