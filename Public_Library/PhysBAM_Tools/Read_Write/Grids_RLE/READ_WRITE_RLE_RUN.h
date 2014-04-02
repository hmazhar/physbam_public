//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RLE_RUN
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RLE_RUN__
#define __READ_WRITE_RLE_RUN__

#include <PhysBAM_Tools/Grids_RLE/RLE_RUN.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<RLE_RUN,RW>
{
public:
    static void Read(std::istream& input,RLE_RUN& object)
    {Read_Binary<RW>(input,object.is_long,object.jmin);}

    static void Write(std::ostream& output,const RLE_RUN& object)
    {Write_Binary<RW>(output,object.is_long,object.jmin);}
};
}

#endif
#endif
#endif
