//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_LEVELSET_1D
//#####################################################################
#ifndef __OPENGL_LEVELSET_1D__
#define __OPENGL_LEVELSET_1D__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_1D.h>
namespace PhysBAM{

template<class T>
class OPENGL_LEVELSET_1D:public OPENGL_SCALAR_FIELD_1D<T,T>
{
public:
    LEVELSET_1D<T>& levelset;
    OPENGL_COLOR point_color, line_color;

    OPENGL_LEVELSET_1D(LEVELSET_1D<T>& levelset_input,const OPENGL_COLOR& point_color_input=OPENGL_COLOR::Cyan(),const OPENGL_COLOR& line_color_input=OPENGL_COLOR::Blue((T).5))
        :OPENGL_SCALAR_FIELD_1D<T,T>(levelset_input.grid,levelset_input.phi,point_color_input,line_color_input),levelset(levelset_input),
        point_color(point_color_input),line_color(line_color_input)
    {}

    ~OPENGL_LEVELSET_1D()
    {}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE
    {
        OPENGL_SCALAR_FIELD_1D<T,T>::Display(in_color);
    }
//#####################################################################
};
}
#endif
