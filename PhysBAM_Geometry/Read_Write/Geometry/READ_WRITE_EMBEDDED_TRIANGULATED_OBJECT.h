//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EMBEDDED_TRIANGULATED_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_EMBEDDED_TRIANGULATED_OBJECT__
#define __READ_WRITE_EMBEDDED_TRIANGULATED_OBJECT__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_EMBEDDED_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<EMBEDDED_TRIANGULATED_OBJECT<TV>,RW>:public Read_Write<EMBEDDED_OBJECT<TV,2>,RW>
{
};
}
#endif
#endif
