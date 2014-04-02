//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RIGID_GEOMETRY_COLLECTION
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RIGID_GEOMETRY_COLLECTION__
#define __READ_WRITE_RIGID_GEOMETRY_COLLECTION__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE_LIST.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,RW>
{
public:
    static void Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame,RIGID_GEOMETRY_COLLECTION<TV>& object,ARRAY<int>* needs_init=0,ARRAY<int>* needs_destroy=0);
    static void Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame,const RIGID_GEOMETRY_COLLECTION<TV>& object);
};
}
#endif
#endif
