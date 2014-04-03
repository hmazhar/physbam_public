//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_INTERVAL
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_INTERVAL__
#define __READ_WRITE_INTERVAL__

#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<INTERVAL<T>,RW>
{
public:
    static void Read(std::istream& input,INTERVAL<T>& object)
    {Read_Binary<RW>(input,object.min_corner,object.max_corner);}

    static void Write(std::ostream& output,const INTERVAL<T>& object)
    {Write_Binary<RW>(output,object.min_corner,object.max_corner);}
};

template<class T>
inline std::ostream& operator<<(std::ostream& output,const INTERVAL<T>& interval)
{output<<"("<<interval.min_corner<<" "<<interval.max_corner<<")";return output;}
}
#endif
#endif
