#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OCTREE_SLICE
//##################################################################### 
#ifndef __OPENGL_OCTREE_SLICE__
#define __OPENGL_OCTREE_SLICE__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE.h>
namespace PhysBAM{

class OPENGL_WORLD;

class OPENGL_OCTREE_SLICE:public OPENGL_SLICE
{
public:
    OPENGL_WORLD &world;
    OCTREE_GRID<float>* octree_grid;

    OPENGL_SLICE::SLICE_MODE mode;

    int axis;
    float position;

    GLenum clip_plane_id1, clip_plane_id2;

    OPENGL_OCTREE_SLICE(OPENGL_WORLD &world_input)
        : world(world_input), octree_grid(0), clip_plane_id1(0), clip_plane_id2(0)
    {}

    bool Is_Slice_Mode() PHYSBAM_OVERRIDE
    {return mode!=NO_SLICE;}

    void Initialize(OCTREE_GRID<float>& octree_grid_input)
    {
        octree_grid=&octree_grid_input;
        Set_Slice_Mode(OPENGL_SLICE::NO_SLICE);
    }

    void Initialize(OCTREE_GRID<double>& octree_grid_input)
    {
        PHYSBAM_FUNCTION_IS_NOT_DEFINED(); // TODO: this class should be templated so this hack can go away...
    }

    void Set_Slice_Mode(OPENGL_SLICE::SLICE_MODE mode_input) PHYSBAM_OVERRIDE
    {
        mode=mode_input;

        if(Is_Slice_Mode()){
            // keep current values of axis and slice number as long as they're in range
            GRID<VECTOR<float,3> >& grid=octree_grid->uniform_grid;
            VECTOR<float,3> lower(grid.domain.min_corner),upper(grid.domain.max_corner);
            float min_value=lower[axis],max_value=upper[axis];
            if(mode==OPENGL_SLICE::CELL_SLICE) max_value-=Get_Slice_Increment();
            position=clamp(position,min_value,max_value);}
        Update_Clip_Planes();
    }

    void Toggle_Slice_Mode() PHYSBAM_OVERRIDE
    {
        axis=clamp(axis,1,3);
        if(mode==OPENGL_SLICE::NO_SLICE){Set_Slice_Mode(OPENGL_SLICE::CELL_SLICE);}
        else if(mode==OPENGL_SLICE::CELL_SLICE){Set_Slice_Mode(OPENGL_SLICE::NODE_SLICE);}
        else if(mode==OPENGL_SLICE::NODE_SLICE){Set_Slice_Mode(OPENGL_SLICE::NO_SLICE);}
    }
    
    void Toggle_Slice_Axis() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()) {
            axis=axis%3+1;
            GRID<VECTOR<float,3> >& grid=octree_grid->uniform_grid;
            position=grid.X(grid.counts/2)[axis];
            Update_Clip_Planes();
        }
    }

    void Increment_Slice() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()) { 
            position+=Get_Slice_Increment();
            float dx=octree_grid->uniform_grid.dX[axis];
            RANGE<VECTOR<float,3> > domain=octree_grid->uniform_grid.domain;
            float max_value=domain.Maximum_Corner()[axis];
            if(mode==OPENGL_SLICE::CELL_SLICE) max_value-=Get_Slice_Increment();
            position=min(position,max_value+octree_grid->number_of_ghost_cells*dx);
            Update_Clip_Planes();
        }
    }

    void Decrement_Slice() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()) {
            position-=Get_Slice_Increment();
            float dx=octree_grid->uniform_grid.dX[axis];
            RANGE<VECTOR<float,3> > domain=octree_grid->uniform_grid.domain;
            position=max(position,domain.Minimum_Corner()[axis]-octree_grid->number_of_ghost_cells*dx);
            Update_Clip_Planes();
        }
    }
    
    void Enable_Clip_Planes() PHYSBAM_OVERRIDE
    {
        glEnable(clip_plane_id1);
        glEnable(clip_plane_id2);
    }

    float Get_Slice_Increment()
    {VECTOR<float,3> DX(octree_grid->uniform_grid.dX);
    return DX[axis]/(1<<(octree_grid->maximum_depth-1));}

    void Print_Slice_Info(std::ostream& output_stream) PHYSBAM_OVERRIDE;
private:
    void Update_Clip_Planes() PHYSBAM_OVERRIDE;
};
}
#endif
#endif
