//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_GREEN_CHILDREN_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_GREEN_CHILDREN_3D__
#define __READ_WRITE_GREEN_CHILDREN_3D__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Red_Green/GREEN_CHILDREN_3D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<GREEN_CHILDREN_3D<T>,RW>
{
public:
    static void Read(std::istream& input,GREEN_CHILDREN_3D<T>& object)
    {int backward_compatible;Read_Binary<RW>(input,backward_compatible,object.elements,object.cells);Read_Binary_Array<RW>(input,object.midpoints,6);} // parent assigned elsewhere

    static void Write(std::ostream& output,const GREEN_CHILDREN_3D<T>& object)
    {Write_Binary<RW>(output,4,object.elements,object.cells);Write_Binary_Array<RW>(output,object.midpoints,6);}
};
}
#endif
#endif
#endif
