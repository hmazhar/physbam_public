//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BINTREE_GRID
//##################################################################### 
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __OPENGL_BINTREE_GRID__
#define __OPENGL_BINTREE_GRID__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>

namespace PhysBAM
{
template<class T>
class OPENGL_BINTREE_GRID:public OPENGL_OBJECT
{
public:
    BINTREE_GRID<T> grid;
    OPENGL_COLOR color;
private:
    OPENGL_SELECTION *current_selection;

public:
    OPENGL_BINTREE_GRID(const OPENGL_COLOR& color_input=OPENGL_COLOR::White()):color(color_input),current_selection(0)
    {}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION *Get_Updated_Selection(OPENGL_SELECTION *selection); // call to update selection (if possible) after grid changes
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_BINTREE_CELL : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_BINTREE_CELL(OPENGL_OBJECT *object,const int index=0)
        :OPENGL_SELECTION(OPENGL_SELECTION::BINTREE_CELL,object),index(index)
    {}
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_BINTREE_NODE : public OPENGL_SELECTION
{
public:
    int index;
    VECTOR<T,1> location;

    OPENGL_SELECTION_BINTREE_NODE(OPENGL_OBJECT *object,const int index,const VECTOR<T,1>& location)
        :OPENGL_SELECTION(OPENGL_SELECTION::BINTREE_NODE,object),index(index),location(location)
    {}
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
#endif
