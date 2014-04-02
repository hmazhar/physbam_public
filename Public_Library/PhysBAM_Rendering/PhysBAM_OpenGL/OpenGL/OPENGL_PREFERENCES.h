//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_PREFERENCES
//#####################################################################
#ifndef __OPENGL_PREFERENCES__
#define __OPENGL_PREFERENCES__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>

namespace PhysBAM{
class OPENGL_PREFERENCES{
public:
    static bool smooth_points;
    static bool smooth_lines;
    static float point_size;
    static float line_width;
    static float selection_point_size;
    static float selection_line_width;
    static float highlighted_point_size;
    static float highlighted_line_width;

    static const float arrowhead_size;
    static const float arrowhead_angle;
    static const OPENGL_COLOR selection_highlight_color;

    static void Set_Smooth_Defaults();
};
}
#endif
