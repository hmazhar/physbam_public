//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TRIANGULATED_AREA
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TRIANGULATED_AREA__
#define __READ_WRITE_TRIANGULATED_AREA__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<TRIANGULATED_AREA<T>,RW>:public Read_Write<MESH_OBJECT<VECTOR<T,2>,TRIANGLE_MESH>,RW>
{
};
}
#endif
#endif
