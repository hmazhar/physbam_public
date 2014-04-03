//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const RANGE<VECTOR<T,2> >& box,const SEGMENT_2D<T>& segment,const T thickness_over_two)
{
    RAY<VECTOR<T,2> > ray(segment);return !box.Outside(segment.x1,thickness_over_two) || !box.Outside(segment.x2,thickness_over_two) || INTERSECTION::Intersects(ray,box,thickness_over_two);
}
//#####################################################################
template bool Intersects(const RANGE<VECTOR<float,2> >&,const SEGMENT_2D<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const RANGE<VECTOR<double,2> >&,const SEGMENT_2D<double>&,const double);
#endif
};
};
