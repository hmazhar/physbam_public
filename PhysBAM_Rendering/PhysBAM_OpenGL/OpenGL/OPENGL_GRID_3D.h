//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Michael Lentine, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_3D
//##################################################################### 
#ifndef __OPENGL_GRID_3D__
#define __OPENGL_GRID_3D__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_GRID_3D : public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    GRID<TV> &grid;
    OPENGL_COLOR color;
    bool draw_ghost_values;
    bool hide_non_selected_grid;
    int scale;

//##################################################################### 
    OPENGL_GRID_3D(GRID<TV> &grid_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White());
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Toggle_Draw_Ghost_Values();
    DEFINE_CALLBACK_CREATOR(OPENGL_GRID_3D,Toggle_Draw_Ghost_Values);
private:
    void Draw_Subgrid(const TV_INT &node_start,const TV_INT &node_end) const;
    void Draw_Nodes_For_Selection(const TV_INT &node_start,const TV_INT &node_end) const;

    OPENGL_SELECTION *current_selection;
//##################################################################### 
};

template<class T>
class OPENGL_SELECTION_GRID_CELL_3D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,3> TV;
public:
    TV_INT index;
    OPENGL_SELECTION_GRID_CELL_3D(OPENGL_OBJECT *object, const TV_INT &index=TV_INT()) 
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_CELL_3D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_GRID_NODE_3D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    TV_INT index;
    OPENGL_SELECTION_GRID_NODE_3D(OPENGL_OBJECT *object, const TV_INT &index=TV_INT())
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_NODE_3D, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_GRID_CELL_LIST_3D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,3> TV;
public:
    ARRAY<TV_INT> indicies;
    OPENGL_SELECTION_GRID_CELL_LIST_3D(OPENGL_OBJECT *object, const ARRAY<TV_INT> &indicies=ARRAY<TV_INT>()) 
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_CELL_LIST_3D, object), indicies(indicies) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_GRID_NODE_LIST_3D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,3> TV;
public:
    ARRAY<TV_INT> indicies;
    OPENGL_SELECTION_GRID_NODE_LIST_3D(OPENGL_OBJECT *object, const ARRAY<TV_INT> &indicies=ARRAY<TV_INT>())
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_NODE_LIST_3D, object), indicies(indicies) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
