//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_MAC_VELOCITY_FIELD_2D
//#####################################################################
#ifndef __OPENGL_MAC_VELOCITY_FIELD_2D__
#define __OPENGL_MAC_VELOCITY_FIELD_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class RANGE;

template<class T_input>
class OPENGL_MAC_VELOCITY_FIELD_2D:public OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef ARRAY<T,FACE_INDEX<2> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    using OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >::size;

    enum VELOCITY_MODE { FACE_CENTERED, CELL_CENTERED };
    VELOCITY_MODE velocity_mode;

    GRID<TV> &grid;
    ARRAY<T,FACE_INDEX<2> > &face_velocities;
    ARRAY_VIEW<T,VECTOR<int,2> > &u,&v;
    ARRAY<TV> vector_field,vector_locations;
    T_ARRAYS_BOOL *active_cells;
    T_FACE_ARRAYS_BOOL *active_faces;

    OPENGL_MAC_VELOCITY_FIELD_2D(GRID<TV> &grid,ARRAY<T,FACE_INDEX<2> > &face_velocities_input,T_ARRAYS_BOOL *active_cells_input=0,T_FACE_ARRAYS_BOOL *active_faces_input=0);
    ~OPENGL_MAC_VELOCITY_FIELD_2D();

    void Update();  // Call when grid/u/v change
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    void Set_Velocity_Mode(VELOCITY_MODE velocity_mode_input);

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    // convenience functions
    void Toggle_Velocity_Mode();
};
}
#endif
