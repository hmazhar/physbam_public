//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D__
#define __OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T,class RW=T>
class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::VECTOR_INT TV_INT;
public:
    OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D(const GRID<TV> &grid,const std::string &velocity_filename_input);
    ~OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_mac_velocity_field.Set_Slice(slice_input);opengl_vorticity_magnitude.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE { opengl_mac_velocity_field.Slice_Has_Changed();opengl_vorticity_magnitude.Slice_Has_Changed(); }
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    void Set_Vector_Size(double size);

    void Toggle_Velocity_Mode();
    void Toggle_Velocity_Mode_And_Draw();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Toggle_Draw_Vorticity();
    void Normalize_Vorticity_Color_Map();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D, Toggle_Velocity_Mode, "Toggle velocity mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D, Toggle_Velocity_Mode_And_Draw, "Toggle velocity mode and draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D, Toggle_Arrowhead, "Toggle arrowhead");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D, Toggle_Draw_Vorticity, "Toggle draw vorticity");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D, Normalize_Vorticity_Color_Map, "Normalize vorticity map based on current frame");

private:
    void Reinitialize();
    void Update_Vorticity();

public:
    OPENGL_MAC_VELOCITY_FIELD_3D<T> opengl_mac_velocity_field;
    OPENGL_SCALAR_FIELD_3D<T> opengl_vorticity_magnitude;
    bool draw_vorticity;

private:
    std::string velocity_filename;
    int frame_loaded;
    bool valid;
    T min_vorticity,max_vorticity;
};

}

#endif
