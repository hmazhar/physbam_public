//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POINT_SIMPLICES_1D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_POINT_SIMPLICES_1D__
#define __READ_WRITE_POINT_SIMPLICES_1D__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
namespace PhysBAM{

void Register_Read_Write_Point_Simplices_1d();

template<class RW,class T>
class Read_Write<POINT_SIMPLICES_1D<T>,RW>:public Read_Write<MESH_OBJECT<VECTOR<T,1>,POINT_SIMPLEX_MESH>,RW>
{
};
}
#endif
#endif
