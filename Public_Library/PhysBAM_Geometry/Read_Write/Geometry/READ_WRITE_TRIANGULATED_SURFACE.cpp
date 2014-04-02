//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TRIANGULATED_SURFACE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
namespace PhysBAM{

void Register_Read_Write_Triangulated_Surface()
{
#define DIMENSION_READ_WRITE_HELPER(T,RW) \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<TRIANGULATED_SURFACE<T> >();

#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
#define RW_HELPER(RW) \
    DIMENSION_READ_WRITE_HELPER(float,RW)
#else
#define RW_HELPER(RW) \
    DIMENSION_READ_WRITE_HELPER(float,RW) \
    DIMENSION_READ_WRITE_HELPER(double,RW)
#endif

    RW_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    RW_HELPER(double);
#endif
}
}
#endif
