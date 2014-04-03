//#####################################################################
// Copyright 2003, 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BOX_3D
//##################################################################### 
#ifndef __OPENGL_BOX_3D__
#define __OPENGL_BOX_3D__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_BOX_3D : public OPENGL_OBJECT
{
public:
    RANGE<VECTOR<T,3> >       &box;
    OPENGL_COLOR    color;

    OPENGL_BOX_3D(RANGE<VECTOR<T,3> > &box_input,const OPENGL_COLOR &color_input = OPENGL_COLOR::White()) 
        : box(box_input), color(color_input)
    {}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
