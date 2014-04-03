//#####################################################################
// Copyright 3009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TWIST
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TWIST__
#define __READ_WRITE_TWIST__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<TWIST<TV>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<TWIST<TV>,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,TWIST<TV>& object)
    {Read_Binary<RW>(input,object.linear,object.angular);}

    static void Write(std::ostream& output,const TWIST<TV>& object)
    {Write_Binary<RW>(output,object.linear,object.angular);}
};
template<class TV> inline std::istream& operator>>(std::istream& input,TWIST<TV>& v)
{FILE_UTILITIES::Ignore(input,'(');input>>v.linear>>v.angular;FILE_UTILITIES::Ignore(input,')');return input;}

template<class TV> inline std::ostream& operator<<(std::ostream& output,const TWIST<TV>& v)
{output<<"("<<v.linear<<"  "<<v.angular<<")";return output;}
}
#endif
#endif
