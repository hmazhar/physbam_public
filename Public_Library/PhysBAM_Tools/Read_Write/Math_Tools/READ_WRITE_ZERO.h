//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ZERO
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ZERO__
#define __READ_WRITE_ZERO__

#include <PhysBAM_Tools/Math_Tools/ZERO.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<ZERO,RW>
{
public:
    static void Read(std::istream& input,ZERO& object)
    {}

    static void Write(std::ostream& output,const ZERO& object)
    {}
};

//#####################################################################
// Stream input and output
//#####################################################################
inline std::ostream&
operator<<(std::ostream& output,const ZERO)
{return output<<0;}
}
#endif
#endif
