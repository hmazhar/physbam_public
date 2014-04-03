//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LOCAL_GRID
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LOCAL_GRID__
#define __READ_WRITE_LOCAL_GRID__

#include <PhysBAM_Tools/Parallel_Computation/LOCAL_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T_GRID>
class Read_Write<LOCAL_GRID<T_GRID>,RW>
{
public:
    static void Read(std::istream& input,LOCAL_GRID<T_GRID>& object)
    {Read_Binary<RW>(input,object.grid);object.Initialize();}
};
}
#endif
#endif
