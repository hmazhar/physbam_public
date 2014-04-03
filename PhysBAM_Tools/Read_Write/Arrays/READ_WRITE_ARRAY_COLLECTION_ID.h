//#####################################################################
// Copyright 2009, Michael Lentine
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ARRAY_COLLECTION_ID
//#####################################################################
#ifndef __READ_WRITE_ARRAY_COLLECTION_ID__
#define __READ_WRITE_ARRAY_COLLECTION_ID__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
namespace PhysBAM{

#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(READ_WRITE_ARRAY_COLLECTION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int READ_WRITE_ARRAY_COLLECTION_ID;
#endif
}
#endif
