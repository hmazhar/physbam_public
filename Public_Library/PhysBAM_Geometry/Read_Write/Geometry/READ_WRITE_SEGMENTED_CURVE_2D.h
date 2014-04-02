//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SEGMENTED_CURVE_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SEGMENTED_CURVE_2D__
#define __READ_WRITE_SEGMENTED_CURVE_2D__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
namespace PhysBAM{

void Register_Read_Write_Segmented_Curve_2d();

template<class RW,class T>
class Read_Write<SEGMENTED_CURVE_2D<T>,RW>:public Read_Write<SEGMENTED_CURVE<VECTOR<T,2> >,RW>
{
};
}
#endif
#endif
