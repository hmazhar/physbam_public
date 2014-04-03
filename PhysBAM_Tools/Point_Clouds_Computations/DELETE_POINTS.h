//#####################################################################
// Copyright 2009, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function DELETE_POINTS
//#####################################################################
#ifndef __DELETE_POINTS__
#define __DELETE_POINTS__

#include <PhysBAM_Tools/Arrays_Computations/RANGE_COMPUTATIONS.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
namespace PhysBAM{

template<class TV> class POINT_CLOUD;
template<class TV> class RANGE;

namespace POINT_CLOUD_COMPUTATIONS{
    template<class TV> int Delete_Points_Outside_Range(POINT_CLOUD<TV>& points,const RANGE<TV>& range)
    {
        ARRAY<int> array;ARRAYS_COMPUTATIONS::Find_Elements_Outside_Range(points.X,range,array);
        int old_count=points.array_collection->Size();
        for(int i=1;i<=array.Size();i++) points.array_collection->Delete_Element(array(i));
        return old_count-points.array_collection->Size();
    }
}
}
#endif
