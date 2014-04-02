//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SPARSE_MATRIX_NXN
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SPARSE_MATRIX_NXN__
#define __READ_WRITE_SPARSE_MATRIX_NXN__

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_NXN.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_SPARSE_VECTOR_ND.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SPARSE_MATRIX_NXN<T>,RW>
{
public:
    static void Read(std::istream& input,SPARSE_MATRIX_NXN<T>& object)
    {object.Clean_Memory();Read_Binary<RW>(input,object.n);object.A.Resize(object.n);
    for(int i=1;i<=object.n;i++){object.A(i)=new SPARSE_VECTOR_ND<T>(object.n);Read_Binary<RW>(input,*object.A(i));}}

    static void Write(std::ostream& output,const SPARSE_MATRIX_NXN<T>& object)
    {Write_Binary<RW>(output,object.n);for(int i=1;i<=object.n;i++)Write_Binary<RW>(output,*object.A(i));}
};
}
#endif
#endif
