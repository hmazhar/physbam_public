//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SYMMETRIC_MATRIX_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_SYMMETRIC_MATRIX_FIELD_2D__
#define __OPENGL_SYMMETRIC_MATRIX_FIELD_2D__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T>
class OPENGL_SYMMETRIC_MATRIX_FIELD_2D:public OPENGL_OBJECT
{
public:
    typedef VECTOR<T,2> TV;GRID<TV> grid;
    const ARRAY<SYMMETRIC_MATRIX<T,2> ,VECTOR<int,2> >& field;
    ARRAY<MATRIX<T,2> ,VECTOR<int,2> > lines;
    ARRAY<PAIR<bool,bool> ,VECTOR<int,2> > positive;
    OPENGL_COLOR positive_color,negative_color;
    T size;

    OPENGL_SYMMETRIC_MATRIX_FIELD_2D(const GRID<TV>& grid_input,const ARRAY<SYMMETRIC_MATRIX<T,2> ,VECTOR<int,2> >& field_input,const T size_input=0.025,
        const OPENGL_COLOR& positive_color_input=OPENGL_COLOR::Red(),const OPENGL_COLOR& negative_color_input=OPENGL_COLOR::Yellow())
        :grid(grid_input),field(field_input),positive_color(positive_color_input),negative_color(negative_color_input),size(size_input)
    {
        if(grid.Is_MAC_Grid())grid=grid.Get_Regular_Grid_At_MAC_Positions();
    }

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual void Update();
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
//##################################################################### 
};
}
#endif
