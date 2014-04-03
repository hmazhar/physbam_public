//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LINE_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LINE_2D__
#define __READ_WRITE_LINE_2D__

#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<LINE_2D<T>,RW>
{
public:
    static void Read(std::istream& input,LINE_2D<T>& object)
    {Read_Binary<RW>(input,object.normal,object.x1);}

    static void Write(std::ostream& output,const LINE_2D<T>& object)
    {Write_Binary<RW>(output,object.normal,object.x1);}
};
}
#endif
#endif
