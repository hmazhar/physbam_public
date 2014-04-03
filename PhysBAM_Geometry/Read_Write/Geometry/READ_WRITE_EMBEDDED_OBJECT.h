//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EMBEDDED_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_EMBEDDED_OBJECT__
#define __READ_WRITE_EMBEDDED_OBJECT__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_OBJECT.h>
namespace PhysBAM{

template<class RW,class TV,int d>
class Read_Write<EMBEDDED_OBJECT<TV,d>,RW>:public Read_Write<STRUCTURE<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object);
    static void Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object);
    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object);
    static void Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object);
};
}
#endif
#endif
