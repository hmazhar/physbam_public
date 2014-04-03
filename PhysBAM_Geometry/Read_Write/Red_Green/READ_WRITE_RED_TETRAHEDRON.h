//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RED_TETRAHEDRON
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RED_TETRAHEDRON__
#define __READ_WRITE_RED_TETRAHEDRON__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Red_Green/RED_TETRAHEDRON.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<RED_TETRAHEDRON<T>,RW>
{
public:
    static void Read(std::istream& input,RED_TETRAHEDRON<T>& object)
    {Read_Binary<RW>(input,object.red_children,object.green_children,object.child_index,object.cell); // owner assigned elsewhere
    if(object.red_children)object.red_children->parent=&object;else if(object.green_children)object.green_children->parent=&object;}

    static void Write(std::ostream& output,const RED_TETRAHEDRON<T>& object)
    {Write_Binary<RW>(output,object.red_children,object.green_children,object.child_index,object.cell);}
};
}
#endif
#endif
#endif
