//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_HEXAHEDRALIZED_VOLUME
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_HEXAHEDRALIZED_VOLUME__
#define __READ_WRITE_HEXAHEDRALIZED_VOLUME__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_HEXAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<HEXAHEDRALIZED_VOLUME<T>,RW>:public Read_Write<MESH_OBJECT<VECTOR<T,3>,HEXAHEDRON_MESH>,RW>
{
};
}
#endif
#endif
