//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SLICE
//#####################################################################
#ifndef __OPENGL_SLICE__
#define __OPENGL_SLICE__

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <iostream>
namespace PhysBAM{

class OPENGL_SLICE
{
public:
    enum SLICE_MODE {NO_SLICE,CELL_SLICE,NODE_SLICE};

    virtual ~OPENGL_SLICE()
    {}

//#####################################################################
    virtual bool Is_Slice_Mode()=0;
    virtual void Set_Slice_Mode(SLICE_MODE mode_input)=0;
    virtual void Toggle_Slice_Mode()=0;
    virtual void Toggle_Slice_Axis()=0;
    virtual void Increment_Slice()=0;
    virtual void Decrement_Slice()=0;
    virtual void Enable_Clip_Planes()=0;
    virtual void Print_Slice_Info(std::ostream& output_stream){}
private:
    virtual void Update_Clip_Planes()=0;
//#####################################################################
};
}
#endif
