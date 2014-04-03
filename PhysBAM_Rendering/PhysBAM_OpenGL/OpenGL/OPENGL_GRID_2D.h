//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_2D
//##################################################################### 
#ifndef __OPENGL_GRID_2D__
#define __OPENGL_GRID_2D__

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
class OPENGL_GRID_2D : public OPENGL_OBJECT
{
public:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef ARRAY<T,FACE_INDEX<2> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    GRID<TV>      &grid;
    T_ARRAYS_BOOL *active_cell_mask,*ghost_cell_mask;
    T_FACE_ARRAYS_BOOL *active_face_mask,*ghost_face_mask;
    T_ARRAYS_BOOL *active_node_mask,*ghost_node_mask;
    OPENGL_COLOR    color;
    bool draw;
    bool draw_ghost_values;
    int draw_mask_type;
private:
    OPENGL_SELECTION *current_selection;
    std::string basedir;
    int frame;

public:
    OPENGL_GRID_2D(GRID<TV> &grid_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White(),const std::string basedir_input="",const int frame_input=0)
        :grid(grid_input),active_cell_mask(0),ghost_cell_mask(0),active_face_mask(0),ghost_face_mask(0),active_node_mask(0),ghost_node_mask(0),color(color_input),draw(true),draw_ghost_values(true),draw_mask_type(0),current_selection(0),basedir(basedir_input),frame(frame_input)
    {Reinitialize();}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual void Set_Frame(int frame_input);
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    void Toggle_Draw_Ghost_Values();
    void Reinitialize();
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    DEFINE_CALLBACK_CREATOR(OPENGL_GRID_2D, Toggle_Draw_Ghost_Values);
};

template<class T>
class OPENGL_SELECTION_GRID_CELL_2D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    VECTOR<int,2> index;
    OPENGL_SELECTION_GRID_CELL_2D(OPENGL_OBJECT *object, const VECTOR<int,2> &index=TV_INT()) 
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_CELL_2D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_GRID_NODE_2D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    VECTOR<int,2> index;
    OPENGL_SELECTION_GRID_NODE_2D(OPENGL_OBJECT *object, const VECTOR<int,2> &index=TV_INT()) 
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_NODE_2D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
