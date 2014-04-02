//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SPARSE_MATRIX_FLAT_NXN
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SPARSE_MATRIX_FLAT_NXN__
#define __READ_WRITE_SPARSE_MATRIX_FLAT_NXN__

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SPARSE_MATRIX_ENTRY<T>,RW>
{
public:
    static void Read(std::istream& input,SPARSE_MATRIX_ENTRY<T>& object)
    {Read_Binary<RW>(input,object.j,object.a);}

    static void Write(std::ostream& output,const SPARSE_MATRIX_ENTRY<T>& object)
    {Write_Binary<RW>(output,object.j,object.a);}
};

template<class RW,class T>
class Read_Write<SPARSE_MATRIX_FLAT_NXN<T>,RW>
{
public:
    static void Read(std::istream& input,SPARSE_MATRIX_FLAT_NXN<T>& object)
    {object.diagonal_index.Clean_Memory();delete object.C;object.C=0;Read_Binary<RW>(input,object.n,object.offsets,object.A);}

    static void Write(std::ostream& output,const SPARSE_MATRIX_FLAT_NXN<T>& object)
    {Write_Binary<RW>(output,object.n,object.offsets,object.A);}
};
}
#endif
#endif
