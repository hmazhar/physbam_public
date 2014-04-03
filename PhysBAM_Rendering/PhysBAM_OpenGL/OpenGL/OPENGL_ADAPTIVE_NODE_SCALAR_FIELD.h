#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD
//##################################################################### 
#ifndef __OPENGL_ADAPTIVE_NODE_SCALAR_FIELD__
#define __OPENGL_ADAPTIVE_NODE_SCALAR_FIELD__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{
template<class T_GRID>
class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD:public OPENGL_OBJECT
{
    typedef typename T_GRID::SCALAR T;
public:
    T_GRID& grid;
    ARRAY<T>& value;
    OPENGL_COLOR_MAP<T>* color_map;
    int point_size;
    bool smooth_shading;

    OPENGL_ADAPTIVE_NODE_SCALAR_FIELD(T_GRID& grid_input,ARRAY<T>& value_input,OPENGL_COLOR_MAP<T>* color_map_input):
        grid(grid_input),value(value_input),color_map(color_map_input),point_size(5),smooth_shading(false)
    {}

    void Set_Smooth_Shading(const bool smooth_shading_input=true)
    {smooth_shading=smooth_shading_input;}

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
