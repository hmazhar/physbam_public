//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_IMPLICIT_OBJECT_TRANSFORMED
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_IMPLICIT_OBJECT_TRANSFORMED.h>
namespace PhysBAM{

void Register_Read_Write_Implicit_Object_Transformed()
{
#define DIMENSION_READ_WRITE_HELPER(T,RW) \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<T,2>,FRAME<VECTOR<T,2> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<T,3>,FRAME<VECTOR<T,3> > > >();

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
