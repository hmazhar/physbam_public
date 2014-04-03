//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D
//#####################################################################
#ifndef __OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D__
#define __OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class T>
class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D:public OPENGL_OBJECT
{
    typedef VECTOR<T,2> TV;
public:
    T height_scale;
    GRID<TV>& grid;
    ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_minus;
    ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_plus;
    LEVELSET_2D<GRID<TV> >& levelset;
    OPENGL_VECTOR_FIELD_3D<T> minus;
    OPENGL_VECTOR_FIELD_3D<T> plus;

    OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D(GRID<TV>& grid,ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_minus,ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_plus,LEVELSET_2D<GRID<TV> >& levelset);
    ~OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D();

//#####################################################################
    void Update();  // Call when grid/V change
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    // convenience functions
    void Scale_Vector_Size(const T scale);
    void Scale_Height(const T height_scale);
//#####################################################################
};
}
#endif
