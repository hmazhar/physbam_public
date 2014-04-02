//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LEVELSET_DYADIC
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LEVELSET_DYADIC__
#define __READ_WRITE_LEVELSET_DYADIC__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_DYADIC.h>
namespace PhysBAM{

template<class RW,class T_GRID>
class Read_Write<LEVELSET_DYADIC<T_GRID>,RW>
{
public:
    static void Read(std::istream& input,LEVELSET_DYADIC<T_GRID>& object)
    {Read_Binary<RW>(input,object.grid,object.phi);
    delete object.normals;object.normals=0;delete object.curvature;object.curvature=0;}

    static void Write(std::ostream& output,const LEVELSET_DYADIC<T_GRID>& object)
    {Write_Binary<RW>(output,object.grid,object.phi);}
};
}
#endif
#endif
#endif
