//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T,RW>& V_minus_component,OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T,RW>& V_plus_component,
       OPENGL_COMPONENT_LEVELSET_2D<T,RW>& levelset_component)
    : OPENGL_COMPONENT("Two-Phase Magnitude Velocity Field 2D"),magnitude_height_scale(0),V_minus_component(V_minus_component),V_plus_component(V_plus_component),levelset_component(levelset_component),
    opengl_two_phase_velocity_magnitude(V_minus_component.opengl_grid_based_vector_field.grid,V_minus_component.opengl_grid_based_vector_field.V,V_plus_component.opengl_grid_based_vector_field.V,levelset_component.opengl_levelset->levelset)
{
    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
~OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D()
{}

template<class T,class RW> bool OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return valid;
}

template<class T,class RW> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Display(const int in_color) const
{
    if (valid && draw) {
        opengl_two_phase_velocity_magnitude.Display(in_color);
    }
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_two_phase_velocity_magnitude.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Toggle_3D_Mode()
{
    if(magnitude_height_scale>0)magnitude_height_scale=0;
    else magnitude_height_scale=(T)1;
    opengl_two_phase_velocity_magnitude.Scale_Height(magnitude_height_scale);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Increase_Point_Size()
{
    opengl_two_phase_velocity_magnitude.Scale_Vector_Size(1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Decrease_Point_Size()
{
    opengl_two_phase_velocity_magnitude.Scale_Vector_Size(1/1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T,RW>::
Reinitialize(const bool force_even_if_not_drawn)
{  
    if(draw||force_even_if_not_drawn){
        V_minus_component.Reinitialize(true);
        V_plus_component.Reinitialize(true);
        levelset_component.Reinitialize(true);
        valid=V_minus_component.Valid_Frame(frame)&&V_plus_component.Valid_Frame(frame)&&levelset_component.Valid_Frame(frame);
        opengl_two_phase_velocity_magnitude.Update();
    }
}

template class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<double,double>;
#endif
