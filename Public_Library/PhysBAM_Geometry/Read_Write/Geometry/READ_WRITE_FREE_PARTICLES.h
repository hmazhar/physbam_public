//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FREE_PARTICLES
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_FREE_PARTICLES__
#define __READ_WRITE_FREE_PARTICLES__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<FREE_PARTICLES<TV>,RW>:public Read_Write<STRUCTURE<TV>,RW>
{
public:
    static void Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {FREE_PARTICLES<TV>& object=dynamic_cast<FREE_PARTICLES<TV>&>(structure_object);
    Read_Binary<RW>(input,object.nodes);}

    static void Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const FREE_PARTICLES<TV>& object=dynamic_cast<const FREE_PARTICLES<TV>&>(structure_object);
    Write_Binary<RW>(output,object.nodes);}
};
}
#endif
#endif
