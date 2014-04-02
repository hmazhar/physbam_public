//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_LEVELSET_3D__
#define __OPENGL_COMPONENT_LEVELSET_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_LEVELSET_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_COMPONENT_LEVELSET_3D(const std::string& levelset_filename,
                                 const std::string& triangulated_surface_filename = "",
                                 const std::string& filename_set_input = "",
                                 const std::string& filename_triangulated_surface_set_input = "",
                                 bool write_generated_triangulated_surface = false,
                                 bool check_triangulated_surface_file_time = true);

    void Set_Surface_Material(const OPENGL_MATERIAL &front_surface_mat,
                              const OPENGL_MATERIAL &back_surface_mat);
    void Set_Overlayed_Surface_Material(const OPENGL_MATERIAL &overlayed_surface_mat);
    void Set_Slice_Color(const OPENGL_COLOR &inside_slice_color,
                         const OPENGL_COLOR &outside_slice_color);

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE;
    virtual void Slice_Has_Changed() { for(int i=1;i<=opengl_levelset_multiviews.m;i++) opengl_levelset_multiviews(i)->Set_Slice(slice); }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    bool Use_Sets() const {return use_sets;}

    void Toggle_Display_Overlay();
    void Toggle_Slice_Color_Mode();
    void Toggle_Smooth_Slice();
    void Next_Set();
    void Previous_Set();
    void Toggle_Draw_Multiple_Levelsets();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_3D, Toggle_Display_Overlay, "Toggle display overlay (in slice mode)");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_3D, Toggle_Slice_Color_Mode, "Toggle solid/gradient slice colors");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_3D, Toggle_Smooth_Slice, "Toggle smooth levelset draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_3D, Next_Set, "Switch to next set");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_3D, Previous_Set, "Switch to previous set");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_3D, Toggle_Draw_Multiple_Levelsets, "Toggle mutliple/single levelset draw");
private:
    void Reinitialize();
    void Reinitialize_Levelset(const std::string& levelset_filename, const std::string& triangulated_surface_filename, OPENGL_LEVELSET_MULTIVIEW<T,RW>* levelset_multiview);

public:
    OPENGL_LEVELSET_MULTIVIEW<T,RW>* opengl_levelset_multiview;
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T,RW>* > opengl_levelset_multiviews;

private:
    std::string levelset_filename;
    std::string triangulated_surface_filename;
    std::string filename_set;
    std::string filename_triangulated_surface_set;
    bool write_generated_triangulated_surface;
    int frame_loaded;
    bool check_triangulated_surface_file_time;
    int set;
    int set_loaded;
    bool use_sets;
    bool draw_multiple_levelsets;
public:
    int ghost_cells;
};

}

#endif
