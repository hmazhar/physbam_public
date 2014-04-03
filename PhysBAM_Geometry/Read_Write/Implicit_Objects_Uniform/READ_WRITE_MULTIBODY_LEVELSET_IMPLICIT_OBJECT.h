//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_MULTIBODY_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_MULTIBODY_LEVELSET_IMPLICIT_OBJECT__
#define __READ_WRITE_MULTIBODY_LEVELSET_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>,RW>:public Read_Write<IMPLICIT_OBJECT<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object) // TODO -- fix to read/write levelsets
    {MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>& object=dynamic_cast<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>&>(structure_object);
    /*Read_Binary(input,levelsets);*/object.Update_Box();object.Update_Minimum_Cell_Size();}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>& object=dynamic_cast<const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>&>(structure_object);
    Write_Binary<RW>(output,*object.levelsets);}
};
}
#endif
#endif
