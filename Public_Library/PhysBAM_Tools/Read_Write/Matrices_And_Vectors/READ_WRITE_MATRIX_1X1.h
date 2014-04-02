//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_MATRIX_1X1
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_MATRIX_1X1__
#define __READ_WRITE_MATRIX_1X1__

#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<MATRIX<T,1>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<MATRIX<T,1>,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,MATRIX<T,1>& object)
    {Read_Binary<RW>(input,object.x11);}

    static void Write(std::ostream& output,const MATRIX<T,1>& object)
    {Write_Binary<RW>(output,object.x11);}
};
}
#endif
#endif
