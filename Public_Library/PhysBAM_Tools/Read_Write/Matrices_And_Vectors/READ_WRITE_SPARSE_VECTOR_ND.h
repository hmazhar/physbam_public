//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SPARSE_VECTOR_ND
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SPARSE_VECTOR_ND__
#define __READ_WRITE_SPARSE_VECTOR_ND__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SPARSE_VECTOR_ND<T>,RW>
{
public:
    static void Read(std::istream& input,SPARSE_VECTOR_ND<T>& object)
    {Read_Binary<RW>(input,object.n,object.number_of_active_indices);
    delete[] object.indices;delete[] object.x;object.indices=new int[object.number_of_active_indices+1];object.x=new T[object.number_of_active_indices+1];
    for(int i=0;i<=object.number_of_active_indices;i++){Read_Binary<RW>(input,object.indices[i]);Read_Binary<RW>(input,object.x[i]);}}

    static void Write(std::ostream& output,const SPARSE_VECTOR_ND<T>& object)
    {Write_Binary<RW>(output,object.n,object.number_of_active_indices);
    for(int i=0;i<=object.number_of_active_indices;i++){Write_Binary<RW>(output,object.indices[i]);Write_Binary<RW>(output,object.x[i]);}}
};
template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const SPARSE_VECTOR_ND<T>& v)
{for(int i=1;i<=v.n;i++)output_stream<<v(i)<<" ";output_stream<<std::endl;return output_stream;}
}
#endif
#endif
