#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RLE_GRID_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RLE_GRID_3D__
#define __OPENGL_COMPONENT_RLE_GRID_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_RLE_GRID_3D:public OPENGL_COMPONENT
{
public:
    OPENGL_RLE_GRID_3D<T> opengl_grid;
private:
    std::string filename;
    int frame_loaded;
    bool valid;
public:

    virtual OPENGL_SELECTION* Get_Selection(GLuint* buffer,int buffer_size)
    {return opengl_grid.Get_Selection(buffer,buffer_size);}

    virtual void Highlight_Selection(OPENGL_SELECTION* selection) PHYSBAM_OVERRIDE
    {opengl_grid.Highlight_Selection(selection);}

    virtual void Clear_Highlight() PHYSBAM_OVERRIDE
    {opengl_grid.Clear_Highlight();}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw && valid;}

//#####################################################################
    OPENGL_COMPONENT_RLE_GRID_3D(const std::string& filename_input);
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
private:
    void Reinitialize();
//#####################################################################
};
}
#endif
#endif
