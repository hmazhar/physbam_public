//#####################################################################
// Copyright 2002, 2003, Robert Bridson, Ronald Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLANE
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
using namespace PhysBAM;
//#####################################################################
// Segment_Plane_Intersection
//#####################################################################
template<class T> bool PLANE<T>::
Segment_Plane_Intersection(const TV& endpoint1,const TV& endpoint2,T& interpolation_fraction) const
{
    T denominator=TV::Dot_Product(endpoint2-endpoint1,normal);
    if(!denominator){interpolation_fraction=FLT_MAX;return false;} // parallel
    interpolation_fraction=TV::Dot_Product(x1-endpoint1,normal)/denominator;
    return (interpolation_fraction>=0 && interpolation_fraction<=1);
}
//#####################################################################
template class PLANE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLANE<double>;
#endif
