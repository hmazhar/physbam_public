//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ARRAY
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ARRAY__
#define __READ_WRITE_ARRAY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_VIEW.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
namespace PhysBAM{

template<class RW,class T,class ID>
class Read_Write<ARRAY<T,ID>,RW>
{
public:
    static void Read(std::istream& input,ARRAY<T,ID>& object)
    {object.Clean_Memory();ID m;Read_Binary<RW>(input,m);
    if(m<ID()) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid negative array size %d",Value(m)));
    if(!m) return;
    object.Exact_Resize(m);
    Read_Binary_Array<RW>(input,object.Get_Array_Pointer(),Value(m));}

    static void Write(std::ostream& output,const ARRAY<T,ID>& object)
    {Write_Prefix(output,object,object.m);}

    static void Write_Prefix(std::ostream& output,const ARRAY<T,ID>& object,const ID prefix)
    {PHYSBAM_ASSERT(ID()<=prefix && prefix<=object.Size());
    Write_Binary<RW>(output,Value(prefix));Write_Binary_Array<RW>(output,object.Get_Array_Pointer(),Value(prefix));}
};
}

#endif
#endif
