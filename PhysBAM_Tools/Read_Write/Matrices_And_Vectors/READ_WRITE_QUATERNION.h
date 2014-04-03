//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_QUATERNION
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_QUATERNION__
#define __READ_WRITE_QUATERNION__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<QUATERNION<T>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<QUATERNION<T>,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,QUATERNION<T>& object)
    {Read_Binary<RW>(input,object.s,object.v);}

    static void Write(std::ostream& output,const QUATERNION<T>& object)
    {Write_Binary<RW>(output,object.s,object.v);}
};
}
#endif
#endif
