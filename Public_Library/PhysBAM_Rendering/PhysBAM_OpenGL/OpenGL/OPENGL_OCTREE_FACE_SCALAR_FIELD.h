#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OCTREE_FACE_SCALAR_FIELD
//##################################################################### 
#ifndef __OPENGL_OCTREE_FACE_SCALAR_FIELD__
#define __OPENGL_OCTREE_FACE_SCALAR_FIELD__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T> class OCTREE_GRID;
template<class T> class OPENGL_COLOR_MAP;

template<class T,class T2=T>
class OPENGL_OCTREE_FACE_SCALAR_FIELD:public OPENGL_OBJECT
{
public:
    OCTREE_GRID<T>& grid;
    ARRAY<T2>& value;
    OPENGL_COLOR_MAP<T2>* color_map;
    int point_size;
    T line_size;
    bool draw_points;

    OPENGL_OCTREE_FACE_SCALAR_FIELD(OCTREE_GRID<T>& grid_input,ARRAY<T2>& value_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool draw_points_input):
        grid(grid_input),value(value_input),color_map(color_map_input),point_size(5),line_size((T).025),draw_points(draw_points_input)
    {}

    void Scale_Vector_Size(const T scale)
    {line_size*=scale;}

//##################################################################### 
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual void Update(){}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
//##################################################################### 
};
}
#endif
#endif
