//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_READ_WRITE_ARRAY_COLLECTION_ID
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_READ_WRITE_ARRAY_COLLECTION_ID__
#define __READ_WRITE_READ_WRITE_ARRAY_COLLECTOIN_ID__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION_ID.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_ELEMENT_ID.h>
namespace PhysBAM{

template<class RW>
class Read_Write<READ_WRITE_ARRAY_COLLECTION_ID,RW>:public Read_Write<ELEMENT_ID<READ_WRITE_ARRAY_COLLECTION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T>,RW>
{
};
}
#endif
#endif
