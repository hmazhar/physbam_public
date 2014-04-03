//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_VECTOR_0D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_VECTOR_0D__
#define __READ_WRITE_VECTOR_0D__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_0D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<VECTOR<T,0>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<VECTOR<T,0>,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,VECTOR<T,0>& object)
    {}

    static void Write(std::ostream& output,const VECTOR<T,0>& object)
    {}
};
//#####################################################################
// Stream input and output
//#####################################################################
template<class T>
inline std::istream& operator>>(std::istream& input,VECTOR<T,0>&)
{if(input.peek()=='[') input.get();if(input.peek()==']') input.get();return input;}
}
#endif
#endif
