//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RED_CHILDREN_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RED_CHILDREN_3D__
#define __READ_WRITE_RED_CHILDREN_3D__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Red_Green/READ_WRITE_RED_TETRAHEDRON.h>
#include <PhysBAM_Geometry/Read_Write/Red_Green/READ_WRITE_TETRAHEDRAL_GROUP.h>
#include <PhysBAM_Geometry/Red_Green/RED_CHILDREN_3D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<RED_CHILDREN_3D<T>,RW>
{
public:
    static void Read(std::istream& input,RED_CHILDREN_3D<T>& object)
    {Read_Binary_Array<RW>(input,object.children,8);for(int i=0;i<8;i++)object.children[i].owner=&object;
    Read_Binary<RW>(input,object.center,object.childrens_size,object.orientation,object.childrens_depth);Read_Binary_Array<RW>(input,object.nodes,10);} // parent assigned elsewhere

    static void Write(std::ostream& output,const RED_CHILDREN_3D<T>& object)
    {Write_Binary_Array<RW>(output,object.children,8);
    Write_Binary<RW>(output,object.center,object.childrens_size,object.orientation,object.childrens_depth);Write_Binary_Array<RW>(output,object.nodes,10);}
};
}
#endif
#endif
#endif
