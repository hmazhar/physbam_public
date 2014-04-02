//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BINTREE_CELL
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_BINTREE_CELL__
#define __READ_WRITE_BINTREE_CELL__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<BINTREE_CELL<T>,RW>
{
public:
    static void Read(std::istream& input,BINTREE_CELL<T>& object);
    static void Read_Helper(std::istream& input,BINTREE_CELL<T>& object);
    static void Write(std::ostream& output,const BINTREE_CELL<T>& object);
    static void Write_Helper(std::ostream& output,const BINTREE_CELL<T>& object);
};
}
#endif
#endif
#endif
