#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_QUADTREE_NODE_VECTOR_FIELD
//##################################################################### 
#ifndef __OPENGL_QUADTREE_NODE_VECTOR_FIELD__
#define __OPENGL_QUADTREE_NODE_VECTOR_FIELD__

#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
namespace PhysBAM{

template<class T_input>
class OPENGL_QUADTREE_NODE_VECTOR_FIELD:public OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    typedef OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > BASE;
    using BASE::size;using BASE::World_Space_Box;

    QUADTREE_GRID<T>& grid;
    ARRAY<VECTOR<T,2> >& V;
    ARRAY<TV> vector_field,vector_locations;

    OPENGL_QUADTREE_NODE_VECTOR_FIELD(QUADTREE_GRID<T>& grid,ARRAY<VECTOR<T,2> >& V)
        :OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >(vector_field,vector_locations),grid(grid),V(V)
    {}

//##################################################################### 
    void Update();
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
//##################################################################### 
};
}
#endif
#endif
