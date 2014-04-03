//#####################################################################
// Copyright 2002, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_RIGID_BODY_HINTS
//##################################################################### 
//
//#####################################################################
// Guendelman - November 2, 2002
//#####################################################################
#ifndef __OPENGL_RIGID_BODY_HINTS__
#define __OPENGL_RIGID_BODY_HINTS__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
namespace PhysBAM
{

class OPENGL_RIGID_BODY_HINTS
{
public:
    OPENGL_MATERIAL material;
    bool            include_bounding_box;

    OPENGL_RIGID_BODY_HINTS()
        : material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR::White())),
          include_bounding_box(true)
    {}

    OPENGL_RIGID_BODY_HINTS(const OPENGL_MATERIAL &material,
                            bool include_bounding_box)
        : material(material), include_bounding_box(include_bounding_box)
    {}

};
}
#include <PhysBAM_Rendering/PhysBAM_OpenGL/Read_Write/OpenGL/READ_WRITE_OPENGL_RIGID_BODY_HINTS.h>
#endif
