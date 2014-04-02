//#####################################################################
// Copyright 2002-2007, Kevin Der, Ronald Fedkiw, Eilene Hao, Sergey Koltakov, Michael Lentine, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SEGMENTED_CURVE_3D
//##################################################################### 
#ifndef __OPENGL_SEGMENTED_CURVE_3D__
#define __OPENGL_SEGMENTED_CURVE_3D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

class OPENGL_SELECTION;

template<class T>
class OPENGL_SEGMENTED_CURVE_3D: public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    const SEGMENTED_CURVE<TV>& curve;
    const OPENGL_SEGMENTED_CURVE_3D<T>* parent_curve;
    mutable ARRAY<int> segment_nodes;
    mutable HASHTABLE<int,TV> vertex_normals;
    OPENGL_COLOR color,color_gray;
    OPENGL_COLOR vertex_color,vertex_position_color;
    bool draw_vertices, draw_vertex_positions, use_solid_color,hide_unselected;

    OPENGL_SEGMENTED_CURVE_3D(const SEGMENTED_CURVE<TV>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan())
        :curve(curve_input),parent_curve(0),color(color_input),color_gray(color_input.Grayscale()),vertex_color(OPENGL_COLOR::Green(0.9)),
        vertex_position_color(OPENGL_COLOR::Magenta()),draw_vertices(false),draw_vertex_positions(false),use_solid_color(true),smooth_normals(false),current_selection(0)
    {}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    OPENGL_SELECTION *Get_Vertex_Selection(int index);
    OPENGL_SELECTION *Get_Segment_Selection(int index);
    OPENGL_SELECTION *Get_Curve_Selection(int index);

    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE;

    void Initialize_Vertex_Normals() const;
    ARRAY<int> Get_Selected_Edges() const;

private:
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;

    bool smooth_normals;

    OPENGL_SELECTION *current_selection;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_3D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_3D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_3D : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_3D(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::SEGMENTED_CURVE_3D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}
#endif
