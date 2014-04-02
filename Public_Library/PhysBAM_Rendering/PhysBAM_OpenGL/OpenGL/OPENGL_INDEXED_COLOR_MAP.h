//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_INDEXED_COLOR_MAP
//#####################################################################
#ifndef __OPENGL_INDEXED_COLOR_MAP__
#define __OPENGL_INDEXED_COLOR_MAP__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>

namespace PhysBAM{

class OPENGL_INDEXED_COLOR_MAP:public OPENGL_COLOR_MAP<int>
{
public:
    enum INDEX_MODE {PERIODIC,CLAMP};
    INDEX_MODE index_mode;
    ARRAY<OPENGL_COLOR,VECTOR<int,1> > color_map;

    OPENGL_INDEXED_COLOR_MAP();
    virtual ~OPENGL_INDEXED_COLOR_MAP();

    void Set_Index_Mode(INDEX_MODE index_mode_input)
    {index_mode=index_mode_input;}

    void Set_Color(int index, const OPENGL_COLOR &color)
    {color_map.Resize(RANGE<VECTOR<int,1> >(min(color_map.domain.min_corner.x,index),max(color_map.domain.max_corner.x,index)));color_map(index)=color;}

    OPENGL_COLOR Lookup(int index) const PHYSBAM_OVERRIDE;
    static OPENGL_INDEXED_COLOR_MAP* Basic_16_Color_Map();
    static OPENGL_INDEXED_COLOR_MAP* Levelset_Multiple_Color_Map();
    static OPENGL_INDEXED_COLOR_MAP* Particle_Multiple_Color_Map();
    static OPENGL_INDEXED_COLOR_MAP* Rigid_Body_Color_Map();
    static OPENGL_INDEXED_COLOR_MAP* Rigid_Body_Back_Color_Map();
};
}
#endif

