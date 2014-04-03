//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_QUADRILATERAL_MESH
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_QUADRILATERAL_MESH__
#define __READ_WRITE_QUADRILATERAL_MESH__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Topology/QUADRILATERAL_MESH.h>
namespace PhysBAM{

template<class RW>
class Read_Write<QUADRILATERAL_MESH,RW>
{
public:
    static void Read(std::istream& input,QUADRILATERAL_MESH& object)
    {Read_Binary<RW>(input,object.number_nodes,object.quadrilaterals);}

    static void Write(std::ostream& output,const QUADRILATERAL_MESH& object)
    {Write_Binary<RW>(output,object.number_nodes,object.quadrilaterals);}
};
}
#endif
#endif
