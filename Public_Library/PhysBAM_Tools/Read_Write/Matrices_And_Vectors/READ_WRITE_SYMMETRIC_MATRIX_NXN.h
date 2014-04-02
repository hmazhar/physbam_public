//#####################################################################
// Copyright 3009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SYMMETRIC_MATRIX_NXN
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SYMMETRIC_MATRIX_NXN__
#define __READ_WRITE_SYMMETRIC_MATRIX_NXN__

#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_NXN.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SYMMETRIC_MATRIX_NXN<T>,RW>
{
public:
    static void Read(std::istream& input,SYMMETRIC_MATRIX_NXN<T>& object)
    {delete[] object.x;
    Read_Binary<RW>(input,object.n);
    assert(object.n>=0);
    object.size=(object.n*object.n+object.n)/2;
    object.x=0;
    if(object.n>0){object.x=new T[object.size];Read_Binary_Array<RW>(input,object.x,object.size);}}

    static void Write(std::ostream& output,const SYMMETRIC_MATRIX_NXN<T>& object)
    {Write_Binary<RW>(output,object.n);Write_Binary_Array<RW>(output,object.x,object.size);}
};
}
#endif
#endif
