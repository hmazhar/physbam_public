//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SIDED_FACE_INDEX
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SIDED_FACE_INDEX__
#define __READ_WRITE_SIDED_FACE_INDEX__

#include <PhysBAM_Tools/Grids_Uniform/SIDED_FACE_INDEX.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
namespace PhysBAM{

template<class RW,int d>
class Read_Write<SIDED_FACE_INDEX<d>,RW>
{
public:
    static void Read(std::istream& input,SIDED_FACE_INDEX<d>& object)
    {Read_Binary<RW>(input,object.side,object.axis,object.index);}

    static void Write(std::ostream& output,const SIDED_FACE_INDEX<d>& object)
    {Write_Binary<RW>(output,object.side,object.axis,object.index);}
};
// global functions
template<int d>
inline std::ostream& operator<<(std::ostream& output,const SIDED_FACE_INDEX<d>& fi)
{output<<"("<<fi.side<<" "<<fi.axis<<" "<<fi.index<<")";return output;}
}
#endif
#endif
