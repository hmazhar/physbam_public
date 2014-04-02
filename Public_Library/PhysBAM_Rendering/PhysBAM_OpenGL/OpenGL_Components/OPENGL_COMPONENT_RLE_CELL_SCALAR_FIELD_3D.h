#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D__
#define __OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_CELL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D:public OPENGL_COMPONENT
{
public:
    OPENGL_RLE_CELL_SCALAR_FIELD_3D<T,T2> opengl_rle_cell_scalar_field;
    mutable OPENGL_TRIANGULATED_SURFACE<T>* opengl_triangulated_surface;
    mutable TRIANGULATED_SURFACE<T>* triangulated_surface;

private:
    std::string scalar_field_filename;
    int frame_loaded;
    bool valid;
    bool draw_surface;
public:
    OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D(RLE_GRID_3D<T> &grid,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map);
    ~OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_rle_cell_scalar_field.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() {if(draw) opengl_rle_cell_scalar_field.Slice_Has_Changed();}
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    void Increase_Point_Size();
    void Decrease_Point_Size();
    void Draw_Surface();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D, Increase_Point_Size, "Increase point size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D, Decrease_Point_Size, "Decrease point size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D, Draw_Surface, "Toggles between slice and surface drawing modes");

    void Reinitialize();
};
}
#endif
#endif
