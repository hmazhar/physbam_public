//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LEVELSET_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LEVELSET_2D__
#define __READ_WRITE_LEVELSET_2D__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class RW,class T_GRID>
class Read_Write<LEVELSET_2D<T_GRID>,RW>:public Read_Write<LEVELSET_UNIFORM<T_GRID>,RW>
{
};
}
#endif
#endif
