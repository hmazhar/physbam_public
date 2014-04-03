//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_CUTTING_SIMPLEX
//#####################################################################
#ifndef __READ_WRITE_CUTTING_SIMPLEX__
#define __READ_WRITE_CUTTING_SIMPLEX__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<CUTTING_SIMPLEX<T,d>,RW>
{
public:
    static void Read(std::istream& input,CUTTING_SIMPLEX<T,d>& object)
    {Read_Binary<RW>(input,object.type,object.parent,object.element_owner,object.nodes,object.weights,object.abs_tol,object.element_original_coordinates,object.simplex_original_coordinates,object.node_in_embedded_simplex);}

    static void Write(std::ostream& output,const CUTTING_SIMPLEX<T,d>& object)
    {Write_Binary<RW>(output,object.type,object.parent,object.element_owner,object.nodes,object.weights,object.abs_tol,object.element_original_coordinates,object.simplex_original_coordinates,object.node_in_embedded_simplex);}
};
}
#endif
