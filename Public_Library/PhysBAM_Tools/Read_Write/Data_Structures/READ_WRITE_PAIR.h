//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_PAIR
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_PAIR__
#define __READ_WRITE_PAIR__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T1,class T2>
class Read_Write<PAIR<T1,T2>,RW>
{
public:

    static void Read(std::istream& input,PAIR<T1,T2>& object)
    {Read_Binary<RW>(input,object.x,object.y);}

    static void Write(std::ostream& output,const PAIR<T1,T2>& object)
    {Write_Binary<RW>(output,object.x,object.y);}
};
template<class T1,class T2>
inline std::istream& operator>>(std::istream& input,PAIR<T1,T2>& p)
{FILE_UTILITIES::Ignore(input,'(');input>>p.x>>p.y;FILE_UTILITIES::Ignore(input,')');return input;}

template<class T1,class T2>
inline std::ostream& operator<<(std::ostream& output,const PAIR<T1,T2>& p)
{output<<"("<<p.x<<" "<<p.y<<")";return output;}
}
#endif
#endif
