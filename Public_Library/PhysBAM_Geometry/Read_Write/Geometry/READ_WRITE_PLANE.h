//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_PLANE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_PLANE__
#define __READ_WRITE_PLANE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<PLANE<T>,RW>
{
public:
    static void Read(std::istream& input,PLANE<T>& object)
    {Read_Binary<RW>(input,object.normal,object.x1);}

    static void Write(std::ostream& output,const PLANE<T>& object)
    {Write_Binary<RW>(output,object.normal,object.x1);}
};
}
#endif
#endif
