//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RLE_GRID_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RLE_GRID_2D__
#define __READ_WRITE_RLE_GRID_2D__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Read_Write/Grids_RLE/READ_WRITE_RLE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<RLE_GRID_2D<T>,RW>:public Read_Write<RLE_GRID<T,RLE_POLICY_2D<T> >,RW>
{
};
}

#endif
#endif
#endif
