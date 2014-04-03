//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POINT_CLOUD_SUBSET
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_POINT_CLOUD_SUBSET__
#define __READ_WRITE_POINT_CLOUD_SUBSET__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_SUBSET.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class TV,class T_POINT_CLOUD>
class Read_Write<POINT_CLOUD_SUBSET<TV,T_POINT_CLOUD>,RW>
{
public:
    static void Read(std::istream& input,POINT_CLOUD_SUBSET<TV,T_POINT_CLOUD>& object)
    {Read_Binary<RW>(input,object.active_indices);object.Update_Subset_Index_From_Element_Index();}

    static void Write(std::ostream& output,const POINT_CLOUD_SUBSET<TV,T_POINT_CLOUD>& object)
    {Write_Binary<RW>(output,object.active_indices);}
};
}
#endif
#endif
