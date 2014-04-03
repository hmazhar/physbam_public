//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_CURVE_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_CURVE_2D__
#define __OPENGL_COMPONENT_CURVE_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_CURVE_2D : public OPENGL_COMPONENT
{
public:
    ARRAY<T,VECTOR<int,1> > x;
    ARRAY<T,VECTOR<int,1> > u;
    ARRAY<T,VECTOR<int,1> > flux;

private:
    std::string x_filename;
    std::string u_filename;
    std::string flux_filename;
    T domain_xmin,domain_xmax;
    T flux_scale;
    bool draw_flux;
    bool draw_du;
    bool draw_piecewise_constant;
    int ghost_nodes;
    int frame_loaded;
    bool valid;
    OPENGL_SELECTION *current_selection;

public:
    OPENGL_COMPONENT_CURVE_2D(const std::string& x_filename_input,
                              const std::string& u_filename_input,
                              const std::string& flux_filename_input,
                              const T domain_xmin_input=0,
                              const T domain_xmax_input=1);
    ~OPENGL_COMPONENT_CURVE_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION *Get_Vertex_Selection(int index);
    OPENGL_SELECTION *Get_Segment_Selection(int index);

    T Area_Under_Curve() const;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Draw_Flux();
    void Toggle_Draw_Du();
    void Toggle_Draw_Piecewise_Constant();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CURVE_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CURVE_2D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CURVE_2D, Toggle_Draw_Flux, "Toggle draw flux");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CURVE_2D, Toggle_Draw_Du, "Toggle draw du");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_CURVE_2D, Toggle_Draw_Piecewise_Constant, "Toggle draw piecewise constant");

private:
    void Reinitialize(bool force=false);
};

template<class T,class RW>
class OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_CURVE_VERTEX_2D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T,class RW>
class OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_CURVE_SEGMENT_2D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
