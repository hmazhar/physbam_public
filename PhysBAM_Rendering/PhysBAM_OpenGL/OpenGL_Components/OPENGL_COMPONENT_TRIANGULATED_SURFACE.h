//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_TRIANGULATED_SURFACE
//##################################################################### 
#ifndef __OPENGL_COMPONENT_TRIANGULATED_SURFACE__
#define __OPENGL_COMPONENT_TRIANGULATED_SURFACE__

#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_TRIANGULATED_SURFACE : public OPENGL_COMPONENT
{
public:
    OPENGL_COMPONENT_TRIANGULATED_SURFACE(const std::string &filename, bool use_display_list = true);
    ~OPENGL_COMPONENT_TRIANGULATED_SURFACE();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size) { return opengl_triangulated_surface.Get_Selection(buffer,buffer_size); }
    virtual void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE { opengl_triangulated_surface.Highlight_Selection(selection); }
    virtual void Clear_Highlight() PHYSBAM_OVERRIDE { opengl_triangulated_surface.Clear_Highlight(); }

private:
    void Reinitialize();    // Needs to be called after some state changes

public:
    TRIANGULATED_SURFACE<T>& triangulated_surface;
    OPENGL_TRIANGULATED_SURFACE<T> opengl_triangulated_surface;

private:
    std::string filename;
    int frame_loaded;
    bool valid;
    bool use_display_list;
};

}

#endif
