//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SCALAR_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_SCALAR_FIELD_3D__
#define __OPENGL_COMPONENT_SCALAR_FIELD_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>

namespace PhysBAM
{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_SCALAR_FIELD_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    // Should be able to combine these two constructors into one (with a default arg) but for some reason I can't get it to compile in linux...
    OPENGL_COMPONENT_SCALAR_FIELD_3D(const GRID<TV> &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input);
    OPENGL_COMPONENT_SCALAR_FIELD_3D(const GRID<TV> &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,
                                     typename OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_MODE draw_mode_input);
    ~OPENGL_COMPONENT_SCALAR_FIELD_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_scalar_field.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() {if(draw) opengl_scalar_field.Slice_Has_Changed();}

    void Toggle_Smooth_Slice();
    void Toggle_Draw_Mode();
    void Toggle_Color_Map();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_3D, Toggle_Smooth_Slice, "Toggle smooth");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_3D, Toggle_Draw_Mode, "Toggle draw mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_3D, Toggle_Color_Map, "Toggle color map");

private:
    void Reinitialize();

public:
    OPENGL_SCALAR_FIELD_3D<T,T2>  opengl_scalar_field;

private:
    std::string scalar_field_filename;
    int frame_loaded;
    bool valid;
};

}

#endif
