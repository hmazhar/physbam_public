//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_UNIFORM_SLICE
//##################################################################### 
#ifndef __OPENGL_UNIFORM_SLICE__
#define __OPENGL_UNIFORM_SLICE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE.h>

namespace PhysBAM{
class OPENGL_WORLD;
class OPENGL_UNIFORM_SLICE:public OPENGL_SLICE
{
public:
    OPENGL_WORLD &world;
    GRID<VECTOR<float,3> > grid;

    OPENGL_SLICE::SLICE_MODE mode;

    int axis;
    int index;

    GLenum clip_plane_id1, clip_plane_id2;

    OPENGL_UNIFORM_SLICE(OPENGL_WORLD &world_input)
        : world(world_input), clip_plane_id1(0), clip_plane_id2(0)
    {
        Initialize(GRID<VECTOR<float,3> >(2,2,2,0,1,0,1,0,1));
    }

    bool Is_Slice_Mode() PHYSBAM_OVERRIDE
    {
        return mode==NODE_SLICE || mode==CELL_SLICE;
    }

    int Maximum_Slice_Index_In_Current_Mode()
    {
        VECTOR<int,3> grid_cells(grid.numbers_of_cells);
        if(mode==NODE_SLICE) return grid_cells[axis]+1;
        else if(mode==CELL_SLICE) return grid_cells[axis];
        else return 0;
    }

    void Initialize(GRID<VECTOR<float,3> > grid_input)
    {
        grid=grid_input;
        mode=NO_SLICE;
        axis=3;
        index=grid_input.counts.z/2;
        Update_Clip_Planes();
    }

    void Set_Slice_Mode(OPENGL_SLICE::SLICE_MODE mode_input) PHYSBAM_OVERRIDE
    {
        mode=mode_input;
        if(Is_Slice_Mode()) {
            // keep current values of axis and slice number as long as they're in range
            axis=clamp(axis,1,3);
            index=clamp(index,1,Maximum_Slice_Index_In_Current_Mode());
        }
        Update_Clip_Planes();
    }

    void Toggle_Slice_Mode() PHYSBAM_OVERRIDE
    {
        Set_Slice_Mode((OPENGL_SLICE::SLICE_MODE)(((int)mode+1)%3));
    }
    
    void Toggle_Slice_Axis() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()) {
            axis=axis%3+1;
            index=Maximum_Slice_Index_In_Current_Mode()/2+1;
            Update_Clip_Planes();
        }
    }

    void Increment_Slice() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()) {
            index=min(index+1,Maximum_Slice_Index_In_Current_Mode());
            Update_Clip_Planes();
        }
    }

    void Decrement_Slice() PHYSBAM_OVERRIDE
    {
        if(Is_Slice_Mode()) {
            index=max(index-1,1);
            Update_Clip_Planes();
        }
    }
    
    void Enable_Clip_Planes() PHYSBAM_OVERRIDE
    {
        glEnable(clip_plane_id1);
        glEnable(clip_plane_id2);
    }

    template<class T>
    static void Get_Face_Index_Range(const OPENGL_UNIFORM_SLICE *slice, const ARRAYS_ND_BASE<VECTOR<T,3> >& array, int face, VECTOR<int,3> &index_start, VECTOR<int,3> &index_end, int scale=1)
    {
        index_start=VECTOR<int,3>(array.domain.min_corner.x,array.domain.min_corner.y,array.domain.min_corner.z);
        index_end=VECTOR<int,3>(array.domain.max_corner.x,array.domain.max_corner.y,array.domain.max_corner.z);
        if(!slice) return;
        else if(slice->mode==CELL_SLICE) {
            if(face==slice->axis) {
                index_start[slice->axis]=max((slice->index-1)/scale+1,index_start[slice->axis]);
                index_end[slice->axis]=min((slice->index-1)/scale+2,index_end[slice->axis]);
            } else {
                index_start[slice->axis]=max((slice->index-1)/scale+1,index_start[slice->axis]);
                index_end[slice->axis]=min((slice->index-1)/scale+1,index_end[slice->axis]);
            }
        } else if(slice->mode==NODE_SLICE) {
            if(face==slice->axis) {
                index_start[slice->axis]=max((slice->index-1)/scale+1,index_start[slice->axis]);
                index_end[slice->axis]=min((slice->index-1)/scale+1,index_end[slice->axis]);
            } else {
                index_start=VECTOR<int,3>(1,1,1);
                index_end=VECTOR<int,3>(0,0,0);
            }
        }
    }

    void Print_Slice_Info(std::ostream& output_stream) PHYSBAM_OVERRIDE;
private:
    void Update_Clip_Planes() PHYSBAM_OVERRIDE;
};
}
#endif
