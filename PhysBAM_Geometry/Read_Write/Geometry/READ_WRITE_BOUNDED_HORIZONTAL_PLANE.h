//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BOUNDED_HORIZONTAL_PLANE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_BOUNDED_HORIZONTAL_PLANE__
#define __READ_WRITE_BOUNDED_HORIZONTAL_PLANE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<BOUNDED_HORIZONTAL_PLANE<TV>,RW>
{
public:
    static void Read(std::istream& input,BOUNDED_HORIZONTAL_PLANE<TV>& object)
    {Read_Binary<RW>(input,object.half_width);}

    static void Write(std::ostream& output,const BOUNDED_HORIZONTAL_PLANE<TV>& object)
    {Write_Binary<RW>(output,object.half_width);}
};
}
#endif
#endif
