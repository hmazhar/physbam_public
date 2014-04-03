//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_INDIRECT_ARRAY
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_INDIRECT_ARRAY__
#define __READ_WRITE_INDIRECT_ARRAY__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T_ARRAY,class T_INDICES>
class Read_Write<INDIRECT_ARRAY<T_ARRAY,T_INDICES>,RW>
{
public:
    static void Read(std::istream& input,INDIRECT_ARRAY<T_ARRAY,T_INDICES>& object)
    {PHYSBAM_NOT_IMPLEMENTED();}

    static void Write(std::ostream& output,const INDIRECT_ARRAY<T_ARRAY,T_INDICES>& object)
    {Write_Binary<RW>(output,object.Size());Write_Binary_Array<RW>(output,object.Get_Array_Pointer(),object.Size());}
};
}
#endif
#endif
