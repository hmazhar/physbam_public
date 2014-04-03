#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_RLE_FACE_SCALAR_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_RLE_FACE_SCALAR_FIELD_3D__
#define __OPENGL_RLE_FACE_SCALAR_FIELD_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_RLE_FACE_SCALAR_FIELD_3D:public OPENGL_OBJECT
{
public:
    typedef RLE_GRID_3D<T> T_GRID;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_X_ITERATOR FACE_X_ITERATOR;
    typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;typedef typename T_GRID::FACE_Z_ITERATOR FACE_Z_ITERATOR;

    RLE_GRID_3D<T>& grid;
    ARRAY<T2>& value;
    OPENGL_COLOR_MAP<T2>* color_map;
    int point_size;
    T line_size;
    bool draw_points;

    mutable T2 max_norm;
    mutable int max_norm_cell;
    mutable VECTOR<int,3> max_norm_I;

    OPENGL_RLE_FACE_SCALAR_FIELD_3D(RLE_GRID_3D<T>& grid_input,ARRAY<T2>& value_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool draw_points_input):
        grid(grid_input),value(value_input),color_map(color_map_input),point_size(5),line_size(0.025),draw_points(draw_points_input)
    {}

    void Scale_Vector_Size(const T scale)
    {line_size*=scale;}

//##################################################################### 
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual void Update(){}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
//##################################################################### 
};
}
#endif
#endif
