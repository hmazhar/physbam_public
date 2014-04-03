//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ARRAY_VIEW
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ARRAY_VIEW__
#define __READ_WRITE_ARRAY_VIEW__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
namespace PhysBAM{

template<class RW,class T,class ID>
class Read_Write<ARRAY_VIEW<T,ID>,RW>
{
public:
    static void Read(std::istream& input,ARRAY_VIEW<T,ID>& object)
    {ID read_size;Read_Binary<RW>(input,read_size);
    if(read_size!=object.Size()) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Expected size %d, read %d",object.Size(),read_size));
    Read_Binary_Array<RW>(input,object.Get_Array_Pointer(),Value(object.Size()));}

    static void Write(std::ostream& output,const ARRAY_VIEW<T,ID>& object)
    {ID m=object.Size();Write_Binary<RW>(output,m);Write_Binary_Array<RW>(output,object.Get_Array_Pointer(),Value(object.m));}
};
}
#endif
#endif
