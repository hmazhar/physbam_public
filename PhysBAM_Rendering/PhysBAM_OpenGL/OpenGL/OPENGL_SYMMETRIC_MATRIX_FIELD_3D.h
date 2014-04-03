//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SYMMETRIC_MATRIX_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_SYMMETRIC_MATRIX_FIELD_3D__
#define __OPENGL_SYMMETRIC_MATRIX_FIELD_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T>
class OPENGL_SYMMETRIC_MATRIX_FIELD_3D:public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    GRID<TV> grid;
    const ARRAY<SYMMETRIC_MATRIX<T,3>,VECTOR<int,3> >& field;
    ARRAY<TRIPLE<VECTOR<T,3>,MATRIX<T,3>,VECTOR<bool,3> > > entries;
    OPENGL_COLOR positive_color,negative_color;
    T size;

    OPENGL_SYMMETRIC_MATRIX_FIELD_3D(const GRID<TV>& grid_input,const ARRAY<SYMMETRIC_MATRIX<T,3>,VECTOR<int,3> >& field_input,const T size_input=0.025,
        const OPENGL_COLOR& positive_color_input=OPENGL_COLOR::Red(),const OPENGL_COLOR& negative_color_input=OPENGL_COLOR::Yellow())
        :grid(grid_input),field(field_input),positive_color(positive_color_input),negative_color(negative_color_input),size(size_input)
    {}

    void Slice_Has_Changed() PHYSBAM_OVERRIDE
    {Update();}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual void Update();
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
//##################################################################### 
};
}
#endif
