//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_GRID_BASED_VECTOR_FIELD_2D__
#define __OPENGL_GRID_BASED_VECTOR_FIELD_2D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
namespace PhysBAM{

template<class T_input>
class OPENGL_GRID_BASED_VECTOR_FIELD_2D:public OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >::size;

    ARRAY<TV> vector_field,vector_locations;
    GRID<TV>& grid;
    ARRAY<VECTOR<T,2>,VECTOR<int,2> >& V;

    OPENGL_GRID_BASED_VECTOR_FIELD_2D(GRID<TV>& grid,ARRAY<VECTOR<T,2>,VECTOR<int,2> >& V);
    ~OPENGL_GRID_BASED_VECTOR_FIELD_2D();

    void Update();  // Call when grid/V change

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

};
}
#endif
