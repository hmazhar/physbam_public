//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FACE_INDEX
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_FACE_INDEX__
#define __READ_WRITE_FACE_INDEX__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
namespace PhysBAM{

template<class RW,int d>
class Read_Write<FACE_INDEX<d>,RW>
{
public:
    static void Read(std::istream& input,FACE_INDEX<d>& object)
    {Read_Binary<RW>(input,object.axis,object.index);}

    static void Write(std::ostream& output,const FACE_INDEX<d>& object)
    {Write_Binary<RW>(output,object.axis,object.index);}
};
// global functions
template<int d>
inline std::ostream& operator<<(std::ostream& output,const FACE_INDEX<d>& fi)
{output<<"("<<fi.axis<<" "<<fi.index<<")";return output;}
}
#endif
#endif
