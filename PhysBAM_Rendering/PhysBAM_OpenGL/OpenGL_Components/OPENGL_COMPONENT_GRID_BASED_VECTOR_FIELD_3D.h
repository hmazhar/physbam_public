//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D__
#define __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_BASED_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T,class RW=T>
class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D(const GRID<TV> &grid,const std::string &vector_field_filename);
    ~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE { slice=slice_input; opengl_grid_based_vector_field.Set_Slice(slice_input); }
    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE { opengl_grid_based_vector_field.Slice_Has_Changed(); }
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Set_Vector_Size(double size);

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D, Toggle_Arrowhead, "Toggle arrowhead");

private:
    void Reinitialize();

public:
    OPENGL_GRID_BASED_VECTOR_FIELD_3D<T> opengl_grid_based_vector_field;

private:
    std::string vector_field_filename;
    int frame_loaded;
    bool valid;
};

}

#endif
