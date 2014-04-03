//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_GREEN_CHILDREN_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_GREEN_CHILDREN_2D__
#define __READ_WRITE_GREEN_CHILDREN_2D__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Red_Green/GREEN_CHILDREN_2D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<GREEN_CHILDREN_2D<T>,RW>
{
public:
    static void Read(std::istream& input,GREEN_CHILDREN_2D<T>& object)
    {int backward_compatible;Read_Binary<RW>(input,backward_compatible,object.elements,object.cells,object.midpoints[0],object.midpoints[1],object.midpoints[2]);} // parent assigned elsewhere

    static void Write(std::ostream& output,const GREEN_CHILDREN_2D<T>& object)
    {Write_Binary<RW>(output,3,object.elements,object.cells,object.midpoints[0],object.midpoints[1],object.midpoints[2]);}
};
}
#endif
#endif
#endif
