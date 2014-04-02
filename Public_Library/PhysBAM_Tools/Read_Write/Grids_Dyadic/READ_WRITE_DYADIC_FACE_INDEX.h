//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_DYADIC_FACE_INDEX
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_DYADIC_FACE_INDEX__
#define __READ_WRITE_DYADIC_FACE_INDEX__

#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_FACE_INDEX.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <iostream>
namespace PhysBAM{

template<class RW,int d>
class Read_Write<DYADIC_FACE_INDEX<d>,RW>
{
public:
    static void Read(std::istream& input,DYADIC_FACE_INDEX<d>& object)
    {Read_Binary<RW>(input,object.axis,object.index,object.first_cell_index,object.second_cell_index);}

    static void Write(std::ostream& output,const DYADIC_FACE_INDEX<d>& object)
    {Write_Binary<RW>(output,object.axis,object.index,object.first_cell_index,object.second_cell_index);}
};
// global functions
template<int d>
inline std::ostream& operator<<(std::ostream& output,const DYADIC_FACE_INDEX<d>& fi)
{output<<"("<<fi.axis<<" "<<fi.index<<")";return output;}
}
#endif
#endif
