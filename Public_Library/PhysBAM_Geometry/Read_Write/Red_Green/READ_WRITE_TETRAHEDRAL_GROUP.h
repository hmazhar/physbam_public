//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TETRAHEDRAL_GROUP
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TETRAHEDRAL_GROUP__
#define __READ_WRITE_TETRAHEDRAL_GROUP__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Red_Green/TETRAHEDRAL_GROUP.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<TETRAHEDRAL_GROUP<T>,RW>
{
public:
    static void Read(std::istream& input,TETRAHEDRAL_GROUP<T>& object)
    {Read_Binary<RW>(input,object.g);}

    static void Write(std::ostream& output,const TETRAHEDRAL_GROUP<T>& object)
    {Write_Binary<RW>(output,object.g);}
};
}
#endif
#endif
#endif
