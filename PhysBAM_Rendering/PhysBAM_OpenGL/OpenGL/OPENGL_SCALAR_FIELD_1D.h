//#####################################################################
// Copyright 2006, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_SCALAR_FIELD_1D__
#define __OPENGL_SCALAR_FIELD_1D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_SCALAR_FIELD_1D:public OPENGL_OBJECT
{
    typedef VECTOR<T,1> TV;
public:
    const GRID<TV>& grid;
    ARRAY<T2,VECTOR<int,1> > &values;
    OPENGL_COLOR point_color,line_color;
    T scale;

    void Scale(const T scale_input)
    {scale=scale*scale_input;}

//#####################################################################
    OPENGL_SCALAR_FIELD_1D(const GRID<TV>& grid,ARRAY<T2,VECTOR<int,1> >& values,OPENGL_COLOR point_color,OPENGL_COLOR line_color);
    virtual ~OPENGL_SCALAR_FIELD_1D();
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
