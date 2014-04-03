//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_FACE_SCALAR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_FACE_SCALAR_FIELD_3D__
#define __OPENGL_FACE_SCALAR_FIELD_3D__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>

namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_FACE_SCALAR_FIELD_3D : public OPENGL_OBJECT
{
public:
    typedef VECTOR<T,3> TV;GRID<TV> grid;
    ARRAY<T2,FACE_INDEX<3> > &face_values;
    ARRAY_VIEW<T2,VECTOR<int,3> > &x_face_values,&y_face_values,&z_face_values;
    OPENGL_COLOR_MAP<T2> *color_map;
    int scale;
private:
    OPENGL_POINTS_3D<T,ARRAY<VECTOR<T,3> > > opengl_points;

//#####################################################################
public:
    OPENGL_FACE_SCALAR_FIELD_3D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<3> > &face_values_input,OPENGL_COLOR_MAP<T2> *color_map_input);
    ~OPENGL_FACE_SCALAR_FIELD_3D();
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Slice_Has_Changed() PHYSBAM_OVERRIDE;    
    virtual void Update();  // Call when values or other attributes have changed
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
