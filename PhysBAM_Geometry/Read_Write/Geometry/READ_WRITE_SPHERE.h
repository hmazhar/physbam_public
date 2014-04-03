//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SPHERE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SPHERE__
#define __READ_WRITE_SPHERE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<SPHERE<TV>,RW>
{
public:
    static void Read(std::istream& input,SPHERE<TV>& object)
    {Read_Binary<RW>(input,object.center,object.radius);}

    static void Write(std::ostream& output,const SPHERE<TV>& object)
    {Write_Binary<RW>(output,object.center,object.radius);}
};
}
#endif
#endif
