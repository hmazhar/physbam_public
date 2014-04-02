//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D__
#define __OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_2D.h>

namespace PhysBAM
{
template<class T,class RW> class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D;
template<class T,class RW> class OPENGL_COMPONENT_LEVELSET_2D;
template<class T,class RW> class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D;

template<class T,class RW=T>
class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D : public OPENGL_COMPONENT
{
public:
    OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T,RW>& V_minus_component,OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T,RW>& V_plus_component,OPENGL_COMPONENT_LEVELSET_2D<T,RW>& levelset_component);
    ~OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Toggle_3D_Mode();
    void Increase_Point_Size();
    void Decrease_Point_Size();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D, Toggle_3D_Mode, "Toggle 3d mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D, Increase_Point_Size, "Increase point size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D, Decrease_Point_Size, "Decrease point size");

private:
    void Reinitialize(const bool force_even_if_not_drawn=false);
    bool valid;

private:
    T magnitude_height_scale;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T,RW>& V_minus_component;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T,RW>& V_plus_component;
    OPENGL_COMPONENT_LEVELSET_2D<T,RW>& levelset_component;
public:
    OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T> opengl_two_phase_velocity_magnitude;
};

}

#endif
