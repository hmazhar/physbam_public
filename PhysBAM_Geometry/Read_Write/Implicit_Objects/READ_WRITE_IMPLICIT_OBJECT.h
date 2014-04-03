//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_IMPLICIT_OBJECT__
#define __READ_WRITE_IMPLICIT_OBJECT__

#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<IMPLICIT_OBJECT<TV>,RW>:public Read_Write<STRUCTURE<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}

    static void Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}

    static void Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}
};
}
#endif
#endif
