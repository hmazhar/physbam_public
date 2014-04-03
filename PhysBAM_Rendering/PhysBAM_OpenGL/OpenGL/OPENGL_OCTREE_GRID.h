#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OCTREE_GRID
//##################################################################### 
#ifndef __OPENGL_OCTREE_GRID__
#define __OPENGL_OCTREE_GRID__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_OCTREE_GRID:public OPENGL_OBJECT
{
public:
    OCTREE_GRID<T> grid;
    OPENGL_COLOR color;
private:
    OPENGL_SELECTION *current_selection;
    int display_list_id;
    mutable bool reinitialize_display_list;
public:

    OPENGL_OCTREE_GRID(const OPENGL_COLOR& color_input=OPENGL_COLOR::White())
        :color(color_input),current_selection(0),display_list_id(-1),reinitialize_display_list(true)
    {display_list_id=glGenLists(1);}

    ~OPENGL_OCTREE_GRID()
    {glDeleteLists(display_list_id,1);}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE {reinitialize_display_list=true;}
    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    void Update(); // Call when grid changes

    OPENGL_SELECTION *Get_Cell_Selection(const int node_index);
    OPENGL_SELECTION *Get_Node_Selection(const int node_index);
    OPENGL_SELECTION *Get_Updated_Selection(OPENGL_SELECTION *selection); // call to update selection (if possible) after grid changes
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;
private:
    void Draw_Coarse_Cells_In_Slice_Mode() const;
    void Draw_Refined_Cells() const;
    
};
 
template<class T>
class OPENGL_SELECTION_OCTREE_CELL : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_OCTREE_CELL(OPENGL_OBJECT *object,const int index=0)
        :OPENGL_SELECTION(OPENGL_SELECTION::OCTREE_CELL,object),index(index)
    {}
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
}; 

template<class T>
class OPENGL_SELECTION_OCTREE_FACE : public OPENGL_SELECTION
{
public:
    int index;
    int node[4];
    int cell[2];
    OPENGL_SELECTION_OCTREE_FACE(OPENGL_OBJECT *object,const int index=0,const int node1=0,const int node2=0,const int node3=0,const int node4=0,const int cell1=0,const int cell2=0)
        :OPENGL_SELECTION(OPENGL_SELECTION::OCTREE_FACE,object),index(index)
    {node[0]=node1;node[1]=node2;node[2]=node3;node[3]=node4;cell[0]=cell1;cell[1]=cell2;}
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
 

template<class T>
class OPENGL_SELECTION_OCTREE_NODE : public OPENGL_SELECTION
{
public:
    int index;
    VECTOR<T,3> location;

    OPENGL_SELECTION_OCTREE_NODE(OPENGL_OBJECT *object,const int index,const VECTOR<T,3>& location)
        :OPENGL_SELECTION(OPENGL_SELECTION::OCTREE_NODE,object),index(index),location(location)
    {}
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};


}

#endif
#endif
