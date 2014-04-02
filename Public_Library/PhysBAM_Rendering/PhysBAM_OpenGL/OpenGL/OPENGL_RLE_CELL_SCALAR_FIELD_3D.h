#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_RLE_CELL_SCALAR_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_RLE_CELL_SCALAR_FIELD_3D__
#define __OPENGL_RLE_CELL_SCALAR_FIELD_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_RLE_CELL_SCALAR_FIELD_3D:public OPENGL_OBJECT
{
public:
    typedef typename RLE_GRID_3D<T>::CELL_ITERATOR CELL_ITERATOR;

    RLE_GRID_3D<T>& grid;
    ARRAY<T2>& value;
    OPENGL_COLOR_MAP<T2>* color_map;
    int point_size;
    bool draw_filled;

    OPENGL_RLE_CELL_SCALAR_FIELD_3D(RLE_GRID_3D<T>& grid_input,ARRAY<T2>& value_input,OPENGL_COLOR_MAP<T2>* color_map_input):
        grid(grid_input),value(value_input),color_map(color_map_input),point_size(5),draw_filled(false)
    {}

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
