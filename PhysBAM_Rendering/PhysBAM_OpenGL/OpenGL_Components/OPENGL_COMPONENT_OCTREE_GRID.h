#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_OCTREE_GRID
//#####################################################################
#ifndef __OPENGL_COMPONENT_OCTREE_GRID__
#define __OPENGL_COMPONENT_OCTREE_GRID__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class T,class RW=T>
class OPENGL_COMPONENT_OCTREE_GRID : public OPENGL_COMPONENT
{
public:
    OPENGL_OCTREE_GRID<T>  opengl_grid;
    std::string filename;
private:
    int frame_loaded;
    bool valid;
public:

    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE {opengl_grid.Slice_Has_Changed();}

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size)
    {return opengl_grid.Get_Selection(buffer,buffer_size);}

    virtual void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE
    {opengl_grid.Highlight_Selection(selection);}

    virtual void Clear_Highlight() PHYSBAM_OVERRIDE
    {opengl_grid.Clear_Highlight();}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw && valid;}

//#####################################################################
    OPENGL_COMPONENT_OCTREE_GRID(const std::string& filename_input);
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;
    void Reinitialize(const bool force=false);
//#####################################################################
};

}

#endif
#endif
