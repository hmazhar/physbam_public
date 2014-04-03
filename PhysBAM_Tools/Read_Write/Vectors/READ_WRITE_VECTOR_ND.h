//#####################################################################
// Copyright 3009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_VECTOR_ND
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_VECTOR_ND__
#define __READ_WRITE_VECTOR_ND__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<VECTOR_ND<T>,RW>
{
public:
    static void Read(std::istream& input,VECTOR_ND<T>& object)
    {assert(object.Owns_Data());Read_Binary<RW>(input,object.n);delete[] object.x;object.x=new T[object.n];Read_Binary_Array<RW>(input,object.x,object.n);}

    static void Write(std::ostream& output,const VECTOR_ND<T>& object)
    {Write_Binary<RW>(output,object.n);Write_Binary_Array<RW>(output,object.x,object.n);}
};
template<class T> inline std::istream& operator>>(std::istream& input,VECTOR_ND<T>& v)
{FILE_UTILITIES::Ignore(input,'[');for(int i=1;i<=v.n;i++) input>>v(i);FILE_UTILITIES::Ignore(input,']');return input;}
}
#endif
#endif
