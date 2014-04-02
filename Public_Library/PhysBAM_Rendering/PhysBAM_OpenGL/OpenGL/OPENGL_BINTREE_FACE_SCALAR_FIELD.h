//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BINTREE_FACE_SCALAR_FIELD
//##################################################################### 
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __OPENGL_BINTREE_FACE_SCALAR_FIELD__
#define __OPENGL_BINTREE_FACE_SCALAR_FIELD__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{
template<class T,class T2=T>
class OPENGL_BINTREE_FACE_SCALAR_FIELD:public OPENGL_OBJECT
{
public:
    BINTREE_GRID<T>& grid;
    ARRAY<T2>& value;
    OPENGL_COLOR_MAP<T2>* color_map;
    int point_size;
    T line_size;
    bool draw_points;

    OPENGL_BINTREE_FACE_SCALAR_FIELD(BINTREE_GRID<T>& grid_input,ARRAY<T2>& value_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool draw_points_input):
        grid(grid_input),value(value_input),color_map(color_map_input),point_size(5),line_size(0.025),draw_points(draw_points_input)
    {}

    void Scale_Vector_Size(const T scale)
    {line_size*=scale;}

    virtual void Update()
    {}

//##################################################################### 
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
//##################################################################### 
};
}
#endif
#endif
