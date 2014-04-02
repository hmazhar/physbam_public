//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ORIENTED_BOX
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ORIENTED_BOX__
#define __READ_WRITE_ORIENTED_BOX__

#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<ORIENTED_BOX<TV>,RW>
{
public:
    static void Read(std::istream& input,ORIENTED_BOX<TV>& object)
    {Read_Binary<RW>(input,object.corner,object.edges);}

    static void Write(std::ostream& output,const ORIENTED_BOX<TV>& object)
    {Write_Binary<RW>(output,object.corner,object.edges);}
};
template<class TV>
inline std::ostream& operator<<(std::ostream& output_stream,const ORIENTED_BOX<TV>& box)
{output_stream<<"("<<box.corner<<")   (";for(int i=1;i<TV::dimension;i++) output_stream<<box.edges.Column(i)<<" : ";output_stream<<box.edges.Column(TV::dimension)<<")";return output_stream;}
}
#endif
#endif
