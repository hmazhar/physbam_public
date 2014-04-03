//#####################################################################
// Copyright 3009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SYMMETRIC_MATRIX_3X3
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SYMMETRIC_MATRIX_3X3__
#define __READ_WRITE_SYMMETRIC_MATRIX_3X3__

#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SYMMETRIC_MATRIX<T,3>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<SYMMETRIC_MATRIX<T,3>,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,SYMMETRIC_MATRIX<T,3>& object)
    {Read_Binary<RW>(input,object.x11,object.x21,object.x31,object.x22,object.x32,object.x33);}

    static void Write(std::ostream& output,const SYMMETRIC_MATRIX<T,3>& object)
    {Write_Binary<RW>(output,object.x11,object.x21,object.x31,object.x22,object.x32,object.x33);}
};
}
#endif
#endif
