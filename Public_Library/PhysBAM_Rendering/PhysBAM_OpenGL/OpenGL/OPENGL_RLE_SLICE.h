#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_RLE_SLICE
//##################################################################### 
#ifndef __OPENGL_RLE_SLICE__
#define __OPENGL_RLE_SLICE__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE.h>

namespace PhysBAM{
class OPENGL_WORLD;
class OPENGL_RLE_SLICE:public OPENGL_SLICE
{
public:
    OPENGL_WORLD &world;
    RLE_GRID_3D<float>* grid;

    OPENGL_SLICE::SLICE_MODE mode;

    int axis;
    int index;

    GLenum clip_plane_id1, clip_plane_id2;

    OPENGL_RLE_SLICE(OPENGL_WORLD &world_input)
        : world(world_input), grid(0), clip_plane_id1(0), clip_plane_id2(0)
    {}

    void Initialize(RLE_GRID_3D<float>& grid_input)
    {
        grid=&grid_input;
        Set_Slice_Mode(OPENGL_SLICE::NO_SLICE);
        axis=0; // special value so we know when we first enter slice mode
    }

    void Initialize(RLE_GRID_3D<double>& grid_input)
    {
        PHYSBAM_FUNCTION_IS_NOT_DEFINED(); // TODO: this class should be templated so this hack can go away...
    }

    bool Is_Slice_Mode() PHYSBAM_OVERRIDE
    {return mode!=NO_SLICE;}

    void Get_Slice_Bounds_In_Current_Mode(int& min_index,int& max_index)
    {
        VECTOR<int,3> grid_cells(grid->uniform_grid.numbers_of_cells+grid->number_of_ghost_cells);
        PHYSBAM_ASSERT(mode==CELL_SLICE);min_index=1-grid->number_of_ghost_cells;max_index=grid_cells[axis];
    }

    void Set_Slice_Mode(OPENGL_SLICE::SLICE_MODE mode_input) PHYSBAM_OVERRIDE
    {
        mode=mode_input;
        if(Is_Slice_Mode()){
            // keep current values of axis and slice number as long as they're in range
            bool first=(axis==0);
            if(first) axis=3; else axis=clamp(axis,1,3);
            int min_index,max_index;Get_Slice_Bounds_In_Current_Mode(min_index,max_index);
            if(first) index=(min_index+max_index)/2; else index=clamp(index,min_index,max_index);}
        Update_Clip_Planes();
    }

    void Toggle_Slice_Mode() PHYSBAM_OVERRIDE
    {
        if(mode==OPENGL_SLICE::NO_SLICE){Set_Slice_Mode(OPENGL_SLICE::CELL_SLICE);}
        else if(mode==OPENGL_SLICE::CELL_SLICE){Set_Slice_Mode(OPENGL_SLICE::NO_SLICE);}
    }
    
    void Toggle_Slice_Axis() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()){
            axis=axis%3+1;
            int min_index,max_index;Get_Slice_Bounds_In_Current_Mode(min_index,max_index);
            index=(min_index+max_index)/2;
            Update_Clip_Planes();}
    }

    void Increment_Slice() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()){ 
            int min_index,max_index;Get_Slice_Bounds_In_Current_Mode(min_index,max_index);
            index=min(index+1,max_index);
            Update_Clip_Planes();}
    }

    void Decrement_Slice() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()){
            int min_index,max_index;Get_Slice_Bounds_In_Current_Mode(min_index,max_index);
            index=max(index-1,min_index);
            Update_Clip_Planes();}
    }
    
    void Enable_Clip_Planes() PHYSBAM_OVERRIDE
    {
        glEnable(clip_plane_id1);
        glEnable(clip_plane_id2);
    }

    void Print_Slice_Info(std::ostream& output_stream) PHYSBAM_OVERRIDE;
private:
    void Update_Clip_Planes() PHYSBAM_OVERRIDE;
};
}
#endif
#endif
