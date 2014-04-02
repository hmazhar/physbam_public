//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RING
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RING__
#define __READ_WRITE_RING__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<RING<T>,RW>
{
public:
    static void Read(std::istream& input,RING<T>& object)
    {Read_Binary<RW>(input,object.plane1.x1,object.plane2.x1,object.outer_radius,object.inner_radius);
    object.Set_Endpoints(object.plane1.x1,object.plane2.x1);}

    static void Write(std::ostream& output,const RING<T>& object)
    {Write_Binary<RW>(output,object.plane1.x1,object.plane2.x1,object.outer_radius,object.inner_radius);}
};
}
#endif
#endif
