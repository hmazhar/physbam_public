//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_DYADIC_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_DYADIC_IMPLICIT_OBJECT__
#define __READ_WRITE_DYADIC_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_OCTREE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_QUADTREE_GRID.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<DYADIC_IMPLICIT_OBJECT<TV>,RW>:public Read_Write<IMPLICIT_OBJECT<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {DYADIC_IMPLICIT_OBJECT<TV>& object=dynamic_cast<DYADIC_IMPLICIT_OBJECT<TV>&>(structure_object);
    Read_Binary<RW>(input,object.levelset);object.Update_Box();object.Update_Minimum_Cell_Size(object.levelset.grid.maximum_depth);object.Update_Phi_Nodes();}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const DYADIC_IMPLICIT_OBJECT<TV>& object=dynamic_cast<const DYADIC_IMPLICIT_OBJECT<TV>&>(structure_object);
    Write_Binary<RW>(output,object.levelset);}
};
}
#endif
#endif
#endif
