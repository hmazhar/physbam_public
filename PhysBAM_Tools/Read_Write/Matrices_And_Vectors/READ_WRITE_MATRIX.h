//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_MATRIX
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_MATRIX__
#define __READ_WRITE_MATRIX__

#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_0X0.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_1X1.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T,int m,int n>
class Read_Write<MATRIX<T,m,n>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<MATRIX<T,m,n>,RW>::value || (m==n && m<=1)>::TYPE>
{
public:
    static void Read(std::istream& input,MATRIX<T,m,n>& object)
    {Read_Binary_Array<RW>(input,object.x,m*n);}

    static void Write(std::ostream& output,const MATRIX<T,m,n>& object)
    {Write_Binary_Array<RW>(output,object.x,m*n);}
};
}
#endif
#endif
