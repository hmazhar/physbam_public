#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_RLE_GRID_3D
//##################################################################### 
#ifndef __OPENGL_RLE_GRID_3D__
#define __OPENGL_RLE_GRID_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_RLE_GRID_3D:public OPENGL_OBJECT
{
public:
    typedef typename RLE_GRID_3D<T>::CELL_ITERATOR CELL_ITERATOR;typedef typename RLE_GRID_3D<T>::FACE_Y_ITERATOR FACE_Y_ITERATOR;

    RLE_GRID_3D<T> grid;
    OPENGL_COLOR color;
private:
    OPENGL_SELECTION* current_selection;
public:

    OPENGL_RLE_GRID_3D(const OPENGL_COLOR& color_input=OPENGL_COLOR::White()):color(color_input),current_selection(0)
    {}

//##################################################################### 
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;
//##################################################################### 
};

template<class T>
class OPENGL_SELECTION_RLE_CELL_3D:public OPENGL_SELECTION
{
public:
    VECTOR<int,3> I;const RLE_RUN_3D* run;
    OPENGL_SELECTION_RLE_CELL_3D(OPENGL_OBJECT* object,const VECTOR<int,3>& I_input,const RLE_RUN_3D* run)
        :OPENGL_SELECTION(OPENGL_SELECTION::RLE_CELL_3D,object),I(I_input),run(run)
    {}
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

//##################################################################### 
}
#endif
#endif
