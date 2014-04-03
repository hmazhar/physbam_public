//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LEVELSET_UNIFORM
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LEVELSET_UNIFORM__
#define __READ_WRITE_LEVELSET_UNIFORM__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class RW,class T_GRID>
class Read_Write<LEVELSET_UNIFORM<T_GRID>,RW>
{
public:
//#####################################################################
    static void Read(std::istream& input,LEVELSET_UNIFORM<T_GRID>& object);
    static void Write(std::ostream& output,const LEVELSET_UNIFORM<T_GRID>& object);
//#####################################################################
};
}
#endif
#endif
