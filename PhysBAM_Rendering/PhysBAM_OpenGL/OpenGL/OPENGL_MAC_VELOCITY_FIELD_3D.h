//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_MAC_VELOCITY_FIELD_3D
//#####################################################################
#ifndef __OPENGL_MAC_VELOCITY_FIELD_3D__
#define __OPENGL_MAC_VELOCITY_FIELD_3D__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class RANGE;

template<class T_input>
class OPENGL_MAC_VELOCITY_FIELD_3D:public OPENGL_VECTOR_FIELD_3D<T_input>
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    using OPENGL_VECTOR_FIELD_3D<T>::slice;
    using OPENGL_VECTOR_FIELD_3D<T>::vector_field;using OPENGL_VECTOR_FIELD_3D<T>::vector_locations;
    using OPENGL_VECTOR_FIELD_3D<T>::size;

    enum VELOCITY_MODE { FACE_CENTERED, CELL_CENTERED };
    VELOCITY_MODE velocity_mode;

    int         max_vectors_3d;
    int         scale;

    GRID<TV> &grid;
    T_FACE_ARRAYS_SCALAR& face_velocities;
    ARRAY_VIEW<T,VECTOR<int,3> > &u, &v, &w;

    OPENGL_MAC_VELOCITY_FIELD_3D(GRID<TV> &grid);
    ~OPENGL_MAC_VELOCITY_FIELD_3D();

    void Update();  // Call when grid/u/v/w change
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    void Set_Velocity_Mode(VELOCITY_MODE velocity_mode_input);

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE { Update(); }

    // convenience functions
    void Toggle_Velocity_Mode();
};
}
#endif
