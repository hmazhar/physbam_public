//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LEVELSET_IMPLICIT_OBJECT__
#define __READ_WRITE_LEVELSET_IMPLICIT_OBJECT__

#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<LEVELSET_IMPLICIT_OBJECT<TV>,RW>:public Read_Write<IMPLICIT_OBJECT<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {LEVELSET_IMPLICIT_OBJECT<TV>& object=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>&>(structure_object);
    Read_Binary<RW>(input,object.levelset);object.Update_Box();object.Update_Minimum_Cell_Size();}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const LEVELSET_IMPLICIT_OBJECT<TV>& object=dynamic_cast<const LEVELSET_IMPLICIT_OBJECT<TV>&>(structure_object);
    Write_Binary<RW>(output,object.levelset);}
};
}
#endif
#endif
