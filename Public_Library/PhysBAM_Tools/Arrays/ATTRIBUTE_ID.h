//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ATTRIBUTE_ID
//#####################################################################
#ifndef __ATTRIBUTE_ID__
#define __ATTRIBUTE_ID__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
namespace PhysBAM{

#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(ATTRIBUTE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(ATTRIBUTE_INDEX,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int ATTRIBUTE_ID;
typedef int ATTRIBUTE_INDEX;
#endif
}
#endif
