//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_GRID
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_GRID__
#define __READ_WRITE_GRID__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<GRID<TV>,RW>
{
public:
    static void Read(std::istream& input,GRID<TV>& object)
    {Read_Binary<RW>(input,object.counts);
    Read_Binary<RW>(input,object.domain,object.MAC_offset);
    object.Initialize(object.counts,object.domain,object.MAC_offset!=0);}

    static void Write(std::ostream& output,const GRID<TV>& object)
    {Write_Binary<RW>(output,object.counts,object.domain,object.MAC_offset);}
};
// global functions
template<class TV>
inline std::ostream& operator<<(std::ostream& output,const GRID<TV>& grid)
{output<<"("<<grid.domain<<" "<<grid.counts<<")";return output;}
}
#endif
#endif
