//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TORUS
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TORUS__
#define __READ_WRITE_TORUS__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<TORUS<T>,RW>
{
public:
    static void Read(std::istream& input,TORUS<T>& object)
    {Read_Binary<RW>(input,object.center,object.axis,object.inner_radius,object.outer_radius);}

    static void Write(std::ostream& output,const TORUS<T>& object)
    {Write_Binary<RW>(output,object.center,object.axis,object.inner_radius,object.outer_radius);}
};
}
#endif
#endif
