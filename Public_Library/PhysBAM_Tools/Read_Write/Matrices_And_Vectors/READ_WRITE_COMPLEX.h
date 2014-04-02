//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_COMPLEX
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_COMPLEX__
#define __READ_WRITE_COMPLEX__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<COMPLEX<T>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<COMPLEX<T>,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,COMPLEX<T>& object)
    {Read_Binary<RW>(input,object.re,object.im);}

    static void Write(std::ostream& output,const COMPLEX<T>& object)
    {Write_Binary<RW>(output,object.re,object.im);}
};
template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const COMPLEX<T>& c)
{output_stream<<"("<<c.re<<" "<<c.im<<")";return output_stream;}
}
#endif
#endif
