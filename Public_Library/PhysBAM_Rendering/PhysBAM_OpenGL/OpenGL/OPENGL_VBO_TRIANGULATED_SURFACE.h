//#####################################################################
// Copyright 2004, Jiayi Chong
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
#ifndef __OPENGL_VBO_TRIANGULATED_SURFACE__
#define __OPENGL_VBO_TRIANGULATED_SURFACE__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
//OpenGL extension Loading Library

const int MAX_VBO_NUM = 20;
//Change this number to suit your graphics card.
//Geforce 4s and 3s should stick to 65535, Geforce FX and Radeon 9500+ can go higher
//Higher MAX_VBO_POINT_NUM means less VBOS => Greater batching and faster rendering but
//going too high will not make your graphics card happy.
const int MAX_VBO_POINT_NUM = 65535; 

namespace PhysBAM
{

template<class T>
class OPENGL_VBO_TRIANGULATED_SURFACE :
    public OPENGL_TRIANGULATED_SURFACE<T>
{
public:

    OPENGL_VBO_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& surface_input,bool smooth_normals_input=true,
        const OPENGL_MATERIAL& material_input=OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Cyan())) ;

    OPENGL_VBO_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& surface_input,bool smooth_normals_input,
        const OPENGL_MATERIAL& front_material_input,const OPENGL_MATERIAL& back_material_input) ;

    ~OPENGL_VBO_TRIANGULATED_SURFACE(void);

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Org_Display(const int in_color=1) const PHYSBAM_OVERRIDE;


    void Create_VBO(void);
    void Use_VBO(void);

protected:
    bool useVBO;
    int useVBOnum;
    int lastVertexfilled;
    int lastNormalfilled;
    GLuint nameVertexVBO[MAX_VBO_NUM];
    GLuint nameNormalVBO[MAX_VBO_NUM];
};

}
#endif
