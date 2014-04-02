//#####################################################################
// Copyright 3009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_VECTOR
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_VECTOR__
#define __READ_WRITE_VECTOR__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_0D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_2D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<VECTOR<T,d>,RW,typename DISABLE_IF<OR<IS_BINARY_IO_SAFE<VECTOR<T,d>,RW>::value,
    OR<d==0,OR<d==1,OR<d==2,d==3>::value>::value>::value>::value>::TYPE>
{
public:
    static void Read(std::istream& input,VECTOR<T,d>& object)
    {Read_Binary_Array<RW>(input,object.array,d);}

    static void Write(std::ostream& output,const VECTOR<T,d>& object)
    {Write_Binary_Array<RW>(output,object.array,d);}
};
//#####################################################################
// Stream input and output
//#####################################################################
template<class T,int d> std::istream& operator>>(std::istream& input,VECTOR<T,d>& v)
{FILE_UTILITIES::Ignore(input,'[');for(int i=0;i<d;i++) input>>v.array[i];FILE_UTILITIES::Ignore(input,']');return input;}
}
#endif
#endif
