//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_1D
//##################################################################### 
#ifndef __OPENGL_GRID_1D__
#define __OPENGL_GRID_1D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>

namespace PhysBAM
{

template<class T>
class OPENGL_GRID_1D : public OPENGL_OBJECT
{
public:
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef ARRAY<T,FACE_INDEX<1> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    GRID<TV> &grid;
    OPENGL_COLOR color;
    bool draw;
    bool draw_ghost_values;
private:
    std::string basedir;
    int frame;

public:
    OPENGL_GRID_1D(GRID<TV> &grid_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White(),const std::string basedir_input="",const int frame_input=0)
        :grid(grid_input),color(color_input),draw(true),draw_ghost_values(true),basedir(basedir_input),frame(frame_input)
    {}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual void Set_Frame(int frame_input);
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Toggle_Draw_Ghost_Values();

    DEFINE_CALLBACK_CREATOR(OPENGL_GRID_1D, Toggle_Draw_Ghost_Values);
};

template<class T>
class OPENGL_SELECTION_GRID_CELL_1D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
public:
    TV_INT index;
    OPENGL_SELECTION_GRID_CELL_1D(OPENGL_OBJECT *object,const TV_INT& index=TV_INT()) 
        :OPENGL_SELECTION(OPENGL_SELECTION::GRID_CELL_1D,object),index(index)
    {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
