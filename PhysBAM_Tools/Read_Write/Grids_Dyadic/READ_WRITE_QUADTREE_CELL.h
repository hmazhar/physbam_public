//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_QUADTREE_CELL
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_QUADTREE_CELL__
#define __READ_WRITE_QUADTREE_CELL__

#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<QUADTREE_CELL<T>,RW>
{
public:
    static void Read(std::istream& input,QUADTREE_CELL<T>& object);
    static void Read_Helper(std::istream& input,QUADTREE_CELL<T>& object);
    static void Write(std::ostream& output,const QUADTREE_CELL<T>& object);
    static void Write_Helper(std::ostream& output,const QUADTREE_CELL<T>& object);
};
}
#endif
#endif
#endif
