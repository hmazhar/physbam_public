#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OCTREE_NODE_VECTOR_FIELD
//##################################################################### 
#ifndef __OPENGL_OCTREE_NODE_VECTOR_FIELD__
#define __OPENGL_OCTREE_NODE_VECTOR_FIELD__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class T_input>
class OPENGL_OCTREE_NODE_VECTOR_FIELD:public OPENGL_VECTOR_FIELD_3D<T_input>
{
    typedef T_input T;
public:
    using OPENGL_VECTOR_FIELD_3D<T>::vector_field;using OPENGL_VECTOR_FIELD_3D<T>::vector_locations;
    using OPENGL_VECTOR_FIELD_3D<T>::size;using OPENGL_VECTOR_FIELD_3D<T>::World_Space_Box;

    OCTREE_GRID<T>& grid;
    ARRAY<VECTOR<T,3> >& V;

    OPENGL_OCTREE_NODE_VECTOR_FIELD(OCTREE_GRID<T>& grid_input,ARRAY<VECTOR<T,3> >& V_input)
        :OPENGL_VECTOR_FIELD_3D<T>(*(new ARRAY<VECTOR<T,3> >),*(new ARRAY<VECTOR<T,3> >)),grid(grid_input),V(V_input)
    {}

//##################################################################### 
    void Slice_Has_Changed() PHYSBAM_OVERRIDE;
    void Update();
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
//##################################################################### 
};
}
#endif
#endif
