//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RLE_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RLE_IMPLICIT_OBJECT__
#define __READ_WRITE_RLE_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Grids_RLE_Level_Sets/READ_WRITE_LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<RLE_IMPLICIT_OBJECT<TV>,RW>:public Read_Write<IMPLICIT_OBJECT<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {RLE_IMPLICIT_OBJECT<TV>& object=dynamic_cast<RLE_IMPLICIT_OBJECT<TV>&>(structure_object);
    Read_Binary<RW>(input,object.levelset);object.Initialize();}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const RLE_IMPLICIT_OBJECT<TV>& object=dynamic_cast<const RLE_IMPLICIT_OBJECT<TV>&>(structure_object);
    Write_Binary<RW>(output,object.levelset);}
};
}
#endif
#endif
#endif
