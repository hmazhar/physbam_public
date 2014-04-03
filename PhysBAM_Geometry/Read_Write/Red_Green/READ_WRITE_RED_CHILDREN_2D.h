//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RED_CHILDREN_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RED_CHILDREN_2D__
#define __READ_WRITE_RED_CHILDREN_2D__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Red_Green/READ_WRITE_RED_TRIANGLE.h>
#include <PhysBAM_Geometry/Red_Green/RED_CHILDREN_2D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<RED_CHILDREN_2D<T>,RW>
{
public:
    static void Read(std::istream& input,RED_CHILDREN_2D<T>& object)
    {Read_Binary_Array<RW>(input,object.children,4);for(int i=0;i<4;i++)object.children[i].owner=&object;
    Read_Binary<RW>(input,object.center,object.childrens_size,object.orientation,object.childrens_depth);Read_Binary_Array<RW>(input,object.nodes,6);} // parent assigned elsewhere

    static void Write(std::ostream& output,const RED_CHILDREN_2D<T>& object)
    {Write_Binary_Array<RW>(output,object.children,4);
    Write_Binary<RW>(output,object.center,object.childrens_size,object.orientation,object.childrens_depth);Write_Binary_Array<RW>(output,object.nodes,6);}
};
}
#endif
#endif
#endif
