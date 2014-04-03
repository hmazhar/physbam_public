//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_BASED_VECTOR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_GRID_BASED_VECTOR_FIELD_3D__
#define __OPENGL_GRID_BASED_VECTOR_FIELD_3D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class TV> class RANGE;

template<class T_input>
class OPENGL_GRID_BASED_VECTOR_FIELD_3D:public OPENGL_VECTOR_FIELD_3D<T_input>
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    using OPENGL_VECTOR_FIELD_3D<T>::slice;
    using OPENGL_VECTOR_FIELD_3D<T>::vector_field;using OPENGL_VECTOR_FIELD_3D<T>::vector_locations;
    using OPENGL_VECTOR_FIELD_3D<T>::size;

    GRID<TV>& grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& V;
    int max_vectors_3d;

    OPENGL_GRID_BASED_VECTOR_FIELD_3D(GRID<TV>& grid,ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& V);
    ~OPENGL_GRID_BASED_VECTOR_FIELD_3D();

    void Update();  // Call when grid/V change

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE { Update(); }
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

};
}
#endif
