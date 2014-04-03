//#####################################################################
// Copyright 2001-2004, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TRIANGULATED_AREA
//##################################################################### 
#ifndef __OPENGL_TRIANGULATED_AREA__
#define __OPENGL_TRIANGULATED_AREA__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_AREA;
template<class T> class OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX;
template<class T> class OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT;
template<class T> class OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE;

template<class T>
class OPENGL_TRIANGULATED_AREA:public OPENGL_OBJECT
{
public:
    TRIANGULATED_AREA<T>& triangulated_area;
    OPENGL_COLOR vertex_color,segment_color,triangle_color,velocity_color;
    OPENGL_SELECTION* current_selection;
    ARRAY<OPENGL_COLOR>* color_map;
    bool draw_vertices,draw_velocities;
    T velocity_scale;

    OPENGL_TRIANGULATED_AREA(TRIANGULATED_AREA<T>& triangulated_area_input,const bool draw_vertices_input=false,
                             const OPENGL_COLOR& vertex_color_input=OPENGL_COLOR::Red(),
                             const OPENGL_COLOR& segment_color_input=OPENGL_COLOR::Green(),
                             const OPENGL_COLOR& triangle_color_input=OPENGL_COLOR::Blue(),
                             ARRAY<OPENGL_COLOR>* color_map_input=0);

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION* Get_Selection(GLuint* buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection,MATRIX<T,3>* transform) const;
    virtual void Set_Color_Map(ARRAY<OPENGL_COLOR>* color_map_input){color_map=color_map_input;}

    OPENGL_SELECTION* Get_Vertex_Selection(int index);
    OPENGL_SELECTION* Get_Segment_Selection(int index);
    OPENGL_SELECTION* Get_Triangle_Selection(int index);

protected:
    void Draw_Vertices() const;
    void Draw_Segments() const;
    void Draw_Triangles(const bool use_color_map=true) const;
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;
    void Draw_Triangles_For_Selection() const;

    friend class OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T>;
    friend class OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T>;
    friend class OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T>;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX:public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX(OPENGL_OBJECT* object,int index=0) 
        :OPENGL_SELECTION(OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX,object),index(index){}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT:public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT(OPENGL_OBJECT* object,int index=0) 
        :OPENGL_SELECTION(OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT,object),index(index){}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE:public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE(OPENGL_OBJECT* object,int index=0) 
        :OPENGL_SELECTION(OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE,object),index(index){}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
