//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTITION_ID
//#####################################################################
#ifndef __PARTITION_ID__
#define __PARTITION_ID__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
namespace PhysBAM{

#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(PARTITION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int PARTITION_ID;
#endif
}
#endif
