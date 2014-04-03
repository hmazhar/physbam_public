//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_POINT_SIMPLEX_1D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const RANGE<VECTOR<T,1> >& box,const POINT_SIMPLEX_1D<T>& point,const T thickness_over_two)
{
    return !box.Outside(point.x1,thickness_over_two);
}
//#####################################################################
// Function Halfspace_Intersection_Size
//#####################################################################
template<class T> T Halfspace_Intersection_Size(const RANGE<VECTOR<T,1> >& box,const POINT_SIMPLEX_1D<T>& halfspace,VECTOR<T,1>* centroid)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
template bool Intersects(const RANGE<VECTOR<float,1> >&,const POINT_SIMPLEX_1D<float>&,const float);
template float Halfspace_Intersection_Size(const RANGE<VECTOR<float,1> >&,const POINT_SIMPLEX_1D<float>&,VECTOR<float,1>*);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const RANGE<VECTOR<double,1> >&,const POINT_SIMPLEX_1D<double>&,const double);
template double Halfspace_Intersection_Size(const RANGE<VECTOR<double,1> >&,const POINT_SIMPLEX_1D<double>&,VECTOR<double,1>*);
#endif
};
};
