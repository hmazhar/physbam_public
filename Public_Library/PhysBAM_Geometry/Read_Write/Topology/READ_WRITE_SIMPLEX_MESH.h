//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SIMPLEX_MESH
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SIMPLEX_MESH__
#define __READ_WRITE_SIMPLEX_MESH__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
namespace PhysBAM{

template<class RW,class T_MESH>
class Read_Write<T_MESH,RW,typename ENABLE_IF<IS_BASE_OF<SIMPLEX_MESH<T_MESH::dimension>,T_MESH>::value>::TYPE>
{
    enum WORKAROUND{d=T_MESH::dimension};
public:
    static void Read(std::istream& input,SIMPLEX_MESH<d>& object);
    static void Write(std::ostream& output,const SIMPLEX_MESH<d>& object);
};
}
#endif
#endif
