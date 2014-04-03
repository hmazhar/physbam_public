//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D__
#define __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_BASED_VECTOR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T,class RW=T>
class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
public:
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D(const GRID<TV> &grid,const std::string &vector_field_filename_input);
    ~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D, Toggle_Arrowhead, "Toggle arrowhead");

private:
    void Reinitialize(bool force_load_even_if_not_drawn=false);
    template<class,class> friend class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D;
public:
    OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>     opengl_grid_based_vector_field;

private:
    std::string vector_field_filename;
    int frame_loaded;
    bool valid;
};

}

#endif
