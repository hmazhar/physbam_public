//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_FACE_SCALAR_FIELD_1D
//#####################################################################
#ifndef __OPENGL_FACE_SCALAR_FIELD_1D__
#define __OPENGL_FACE_SCALAR_FIELD_1D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_2D.h>

namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_FACE_SCALAR_FIELD_1D : public OPENGL_OBJECT
{
    typedef VECTOR<T,1> TV;
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
public:
    GRID<TV> grid;
    ARRAY<T2,FACE_INDEX<1> > &face_values;
    ARRAY_VIEW<T2,VECTOR<int,1> > &x_face_values;
    OPENGL_COLOR point_color;
    OPENGL_COLOR line_color;
private:
    T scale;
    
public:
    void Scale(const T scale_input)
    {scale=scale*scale_input;}

//#####################################################################
    OPENGL_FACE_SCALAR_FIELD_1D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<1> > &face_values_input,OPENGL_COLOR point_color_input,OPENGL_COLOR line_color_input);
    ~OPENGL_FACE_SCALAR_FIELD_1D();
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
