//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D__
#define __OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
public:
    GRID<TV> grid,mac_grid,u_grid,v_grid;
    ARRAY<VECTOR<bool,2> ,VECTOR<int,2> > node_neighbors_visible;
    ARRAY<VECTOR<PAIR<bool,T>,2>,VECTOR<int,2> > face_corners_visible_from_face_center_u; // length 2, order is bottom, top
    ARRAY<VECTOR<PAIR<bool,T>,2>,VECTOR<int,2> > face_corners_visible_from_face_center_v; // length 2, order is left, right
    ARRAY<bool,VECTOR<int,2> > density_valid_mask;
    ARRAY<bool,VECTOR<int,2> > phi_valid_mask;
private:
    OPENGL_CONSTANT_COLOR_MAP<bool> invalid_color_map;
public:
    OPENGL_SCALAR_FIELD_2D<T,bool> opengl_density_valid_mask;
    OPENGL_SCALAR_FIELD_2D<T,bool> opengl_phi_valid_mask;
private:
    std::string directory;
    int frame_loaded;
    bool valid;
    bool draw_grid_visibility,draw_density_valid_mask,draw_phi_valid_mask;

public:
    OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D(GRID<TV> &grid,const std::string& directory);
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Toggle_Draw_Grid_Visibility();
    void Toggle_Draw_Density_Valid_Mask();
    void Toggle_Draw_Phi_Valid_Mask();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D, Toggle_Draw_Grid_Visibility, "Toggle draw grid visibility");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D, Toggle_Draw_Density_Valid_Mask, "Toggle draw density valid mask");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D, Toggle_Draw_Phi_Valid_Mask, "Toggle draw phi valid mask");

private:
    void Reinitialize(bool force=false);
};
}
#endif
