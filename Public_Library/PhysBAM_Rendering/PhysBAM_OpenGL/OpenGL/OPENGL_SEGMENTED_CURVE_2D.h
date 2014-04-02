//#####################################################################
// Copyright 2002, Ronald Fedkiw, Eilene Hao, Sergey Koltakov, Neil Molino, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SEGMENTED_CURVE_2D
//##################################################################### 
#ifndef __OPENGL_SEGMENTED_CURVE_2D__
#define __OPENGL_SEGMENTED_CURVE_2D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

class OPENGL_SELECTION;

template<class T>
class OPENGL_SEGMENTED_CURVE_2D: public OPENGL_OBJECT
{
public:
    const SEGMENTED_CURVE_2D<T>& curve;
    OPENGL_COLOR color,color_gray;
    OPENGL_COLOR vertex_color,vertex_position_color,velocity_color;
    bool draw_vertices,draw_vertex_positions,draw_velocities;
    T velocity_scale;
private:
    OPENGL_SELECTION *current_selection;

public:
    OPENGL_SEGMENTED_CURVE_2D(const SEGMENTED_CURVE_2D<T>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection,MATRIX<T,3>* transform) const;

    OPENGL_SELECTION *Get_Vertex_Selection(int index);
    OPENGL_SELECTION *Get_Segment_Selection(int index);

private:
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}
#endif
