//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SPARSE_MATRIX_FLAT_MXN
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SPARSE_MATRIX_FLAT_MXN__
#define __READ_WRITE_SPARSE_MATRIX_FLAT_MXN__

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SPARSE_MATRIX_FLAT_MXN<T>,RW>
{
public:
    static void Read(std::istream& input,SPARSE_MATRIX_FLAT_MXN<T>& object)
    {Read_Binary<RW>(input,object.m,object.n,object.offsets,object.A);}

    static void Write(std::ostream& output,const SPARSE_MATRIX_FLAT_MXN<T>& object)
    {Write_Binary<RW>(output,object.m,object.n,object.offsets,object.A);}
};
}
#endif
#endif
