//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_BODY_ID
//#####################################################################
#ifndef __COMBINED_BODY_ID__
#define __COMBINED_BODY_ID__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
namespace PhysBAM{

#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(COMBINED_BODY_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int COMBINED_BODY_ID;
#endif
}
#endif
