//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FAST_LEVELSET
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_FAST_LEVELSET__
#define __READ_WRITE_FAST_LEVELSET__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
namespace PhysBAM{

template<class RW,class T_GRID>
class Read_Write<FAST_LEVELSET<T_GRID>,RW>:public Read_Write<typename LEVELSET_POLICY<T_GRID>::LEVELSET,RW>
{
};
}
#endif
#endif
