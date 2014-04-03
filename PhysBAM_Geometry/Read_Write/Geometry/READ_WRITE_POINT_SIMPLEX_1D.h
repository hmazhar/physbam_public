//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POINT_SIMPLEX_1D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_POINT_SIMPLEX_1D__
#define __READ_WRITE_POINT_SIMPLEX_1D__

#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<POINT_SIMPLEX_1D<T>,RW>
{
public:
    static void Read(std::istream& input,POINT_SIMPLEX_1D<T>& object)
    {Read_Binary<RW>(input,object.direction,object.x1);}

    static void Write(std::ostream& output,const POINT_SIMPLEX_1D<T>& object)
    {Write_Binary<RW>(output,object.direction,object.x1);}
};
}
#endif
#endif
