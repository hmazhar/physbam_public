//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POLYGON_MESH
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_POLYGON_MESH__
#define __READ_WRITE_POLYGON_MESH__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Topology/POLYGON_MESH.h>
namespace PhysBAM{

template<class RW>
class Read_Write<POLYGON_MESH,RW>
{
public:
    static void Read(std::istream& input,POLYGON_MESH& object)
    {object.Clean_Memory();Read_Binary<RW>(input,object.number_nodes,object.elements);}

    static void Write(std::ostream& output,const POLYGON_MESH& object)
    {Write_Binary<RW>(output,object.number_nodes,object.elements);}
};
}
#endif
#endif
