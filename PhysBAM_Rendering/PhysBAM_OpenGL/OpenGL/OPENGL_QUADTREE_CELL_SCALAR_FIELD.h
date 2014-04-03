#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2009, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_QUADTREE_CELL_SCALAR_FIELD
//##################################################################### 
#ifndef __OPENGL_QUADTREE_CELL_SCALAR_FIELD__
#define __OPENGL_QUADTREE_CELL_SCALAR_FIELD__

#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_QUADTREE_CELL_SCALAR_FIELD:public OPENGL_OBJECT
{
public:
    QUADTREE_GRID<T>& grid;
    ARRAY<T2>& value;
    OPENGL_COLOR_MAP<T2>* color_map;
    int point_size;

    OPENGL_QUADTREE_CELL_SCALAR_FIELD(QUADTREE_GRID<T>& grid_input,ARRAY<T2>& value_input,OPENGL_COLOR_MAP<T2>* color_map_input):
        grid(grid_input),value(value_input),color_map(color_map_input),point_size(5)
    {}

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
