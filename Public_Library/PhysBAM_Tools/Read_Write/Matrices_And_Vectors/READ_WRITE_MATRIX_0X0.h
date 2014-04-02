//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_MATRIX_0X0
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_MATRIX_0X0__
#define __READ_WRITE_MATRIX_0X0__

#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<MATRIX<T,0>,RW>
{
public:
    static void Read(std::istream& input,MATRIX<T,0>& object)
    {}

    static void Write(std::ostream& output,const MATRIX<T,0>& object)
    {}
};
}
#endif
#endif
