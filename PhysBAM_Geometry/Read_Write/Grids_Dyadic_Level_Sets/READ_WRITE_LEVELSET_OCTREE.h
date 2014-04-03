//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LEVELSET_OCTREE
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LEVELSET_OCTREE__
#define __READ_WRITE_LEVELSET_OCTREE__

#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_OCTREE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_DYADIC.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<LEVELSET_OCTREE<T>,RW>:public Read_Write<LEVELSET_DYADIC<OCTREE_GRID<T> >,RW>
{
};
}
#endif
#endif
#endif
