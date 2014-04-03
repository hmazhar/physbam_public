//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_HEXAHEDRON_MESH
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_HEXAHEDRON_MESH__
#define __READ_WRITE_HEXAHEDRON_MESH__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Topology/HEXAHEDRON_MESH.h>
namespace PhysBAM{

template<class RW>
class Read_Write<HEXAHEDRON_MESH,RW>
{
public:
    static void Read(std::istream& input,HEXAHEDRON_MESH& object)
    {object.Clean_Memory();int backward_compatible;Read_Binary<RW>(input,object.number_nodes,backward_compatible);Read_Binary<RW>(input,object.elements);}

    static void Write(std::ostream& output,const HEXAHEDRON_MESH& object)
    {Write_Binary<RW>(output,object.number_nodes,8);Write_Binary<RW>(output,object.elements);}
};
}
#endif
#endif
