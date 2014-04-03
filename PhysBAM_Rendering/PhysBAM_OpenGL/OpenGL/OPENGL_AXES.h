//#####################################################################
// Copyright 2002-2005, Robert Bridson, Eilene Hao, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_AXES
//##################################################################### 
#ifndef __OPENGL_AXES__
#define __OPENGL_AXES__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
namespace PhysBAM{

template<class T>
class OPENGL_AXES:public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
    typedef FRAME<TV> T_FRAME;
public:
    RANGE<VECTOR<T,3> > box; // extents of axes with respect to local frame
    bool draw_box; // whether to draw a bounding box or axis vectors
    bool draw_xz_grid,draw_xy_grid,draw_yz_grid; // whether to draw grids on each plane
    T grid_spacing;

    OPENGL_AXES(const FRAME<VECTOR<T,3> >& frame_input=T_FRAME(),const RANGE<VECTOR<T,3> >& box_input=(RANGE<VECTOR<T,3> >(0,1,0,1,0,1)),
        bool draw_box_input=false,bool draw_xz_grid_input=false,bool draw_xy_grid_input=false,bool draw_yz_grid_input=false,T grid_spacing_input=.1)
        :box(box_input),draw_box(draw_box_input),draw_xz_grid(draw_xz_grid_input),draw_xy_grid(draw_xy_grid_input),draw_yz_grid(draw_yz_grid_input),grid_spacing(grid_spacing_input)
    {
        *frame=FRAME<VECTOR<float,3> >(frame_input);
    }

    void Scale(const T scale)
    {box*=scale;}

//#####################################################################
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    bool Is_Transparent() const PHYSBAM_OVERRIDE {return false;}
//#####################################################################
};
}
#endif
