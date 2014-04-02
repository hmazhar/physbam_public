//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BOX
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_BOX__
#define __READ_WRITE_BOX__

#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<BOX<TV>,RW>:public Read_Write<RANGE<TV>,RW>
{
};
}
#endif
#endif
