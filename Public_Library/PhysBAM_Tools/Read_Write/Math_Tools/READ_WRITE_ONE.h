//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ONE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ONE__
#define __READ_WRITE_ONE__

#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<ONE,RW>
{
public:
    static void Read(std::istream& input,ONE& object)
    {}

    static void Write(std::ostream& output,const ONE& object)
    {}
};

//#####################################################################
// Stream input and output
//#####################################################################
inline std::ostream&
operator<<(std::ostream& output,const ONE)
{return output<<1;}
}
#endif
#endif
