//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RLE_RUN_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RLE_RUN_3D__
#define __READ_WRITE_RLE_RUN_3D__

#include <PhysBAM_Tools/Grids_RLE/RLE_RUN_3D.h>
#include <PhysBAM_Tools/Read_Write/Grids_RLE/READ_WRITE_RLE_RUN.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<RLE_RUN_3D,RW>:public Read_Write<RLE_RUN,RW>
{
};
}

#endif
#endif
#endif
