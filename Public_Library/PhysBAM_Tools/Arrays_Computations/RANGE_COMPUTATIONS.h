//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE_COMPUTATIONS
//#####################################################################
#ifndef __RANGE_COMPUTATIONS__
#define __RANGE_COMPUTATIONS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

namespace ARRAYS_COMPUTATIONS
{
    template<class T_ARRAY> static RANGE<typename T_ARRAY::ELEMENT> Compute_Range(const T_ARRAY& array)
    {typedef typename T_ARRAY::ELEMENT T;
    RANGE<T> range;range.min_corner=array(1);range.max_corner=array(1);
    for(int i=2;i<=array.Size();i++){if(array(i)<range.min_corner) range.min_corner=array(i);if(array(i)>range.max_corner) range.max_corner=array(i);}
    return range;}

    template<class T_ARRAY> static void Find_Elements_Inside_Range(const T_ARRAY& array,const RANGE<typename T_ARRAY::ELEMENT>& range,ARRAY<typename T_ARRAY::INDEX>& inside)
    {for(int i=1;i<=array.Size();i++) if(range.Lazy_Inside(array(i))) inside.Append(i);}

    template<class T_ARRAY> static void Find_Elements_Outside_Range(const T_ARRAY& array,const RANGE<typename T_ARRAY::ELEMENT>& range,ARRAY<typename T_ARRAY::INDEX>& outside)
    {for(int i=1;i<=array.Size();i++) if(range.Lazy_Inside(array(i))) outside.Append(i);}
}
}
#endif
