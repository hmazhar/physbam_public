//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TRIPLE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TRIPLE__
#define __READ_WRITE_TRIPLE__

#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T1,class T2,class T3>
class Read_Write<TRIPLE<T1,T2,T3>,RW>
{
public:

    static void Read(std::istream& input,TRIPLE<T1,T2,T3>& object)
    {Read_Binary<RW>(input,object.x,object.y,object.z);}

    static void Write(std::ostream& output,const TRIPLE<T1,T2,T3>& object)
    {Write_Binary<RW>(output,object.x,object.y,object.z);}
};

template<class T1,class T2,class T3>
inline std::ostream& operator<<(std::ostream& output,const TRIPLE<T1,T2,T3>& object)
{output<<"("<<object.x<<" "<<object.y<<" "<<object.z<<")";return output;}
}
#endif
#endif
