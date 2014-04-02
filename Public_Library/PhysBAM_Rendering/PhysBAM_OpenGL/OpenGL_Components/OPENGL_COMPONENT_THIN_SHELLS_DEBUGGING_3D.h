//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D__
#define __OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    GRID<TV> grid,mac_grid,u_grid,v_grid,w_grid;
    ARRAY<VECTOR<bool,3>,VECTOR<int,3> > node_neighbors_visible;
    ARRAY<VECTOR<PAIR<bool,T>,4>,VECTOR<int,3> > face_corners_visible_from_face_center_u; // length 4, order is front bottom, front top, back bottom, back top
    ARRAY<VECTOR<PAIR<bool,T>,4>,VECTOR<int,3> > face_corners_visible_from_face_center_v; // length 4, order is front left, front right, back left, back right
    ARRAY<VECTOR<PAIR<bool,T>,4>,VECTOR<int,3> > face_corners_visible_from_face_center_w; // length 4, order is bottom left, bottom right, top left, top right
    ARRAY<bool,VECTOR<int,3> > density_valid_mask;
private:
    OPENGL_CONSTANT_COLOR_MAP<bool> invalid_color_map;
public:
    OPENGL_SCALAR_FIELD_3D<T,bool> opengl_density_valid_mask;
private:
    std::string directory;
    int frame_loaded;
    bool valid;
    bool draw_density_valid_mask;
    bool draw_node_neighbors_visible;
    bool draw_face_corners_visible;

public:
    OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D(const GRID<TV> &grid,const std::string& directory);
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE;
    void Slice_Has_Changed() PHYSBAM_OVERRIDE;    

    void Toggle_Draw_Grid_Visibility_Mode();
    void Toggle_Draw_Density_Valid_Mask();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D, Toggle_Draw_Grid_Visibility_Mode, "Toggle draw grid visibility mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D, Toggle_Draw_Density_Valid_Mask, "Toggle draw density valid mask");

private:
    void Reinitialize(bool force=false);
};
}
#endif
