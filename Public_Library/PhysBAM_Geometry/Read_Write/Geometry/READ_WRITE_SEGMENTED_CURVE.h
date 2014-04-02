//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SEGMENTED_CURVE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SEGMENTED_CURVE__
#define __READ_WRITE_SEGMENTED_CURVE__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<SEGMENTED_CURVE<TV>,RW>:public Read_Write<MESH_OBJECT<TV,SEGMENT_MESH>,RW>
{
};
}
#endif
#endif
