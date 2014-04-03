//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_POINT_SIMPLEX_1D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLEX_1D<T>& point,const T thickness_over_two)
{
    return Intersects(ray,RANGE<VECTOR<T,1> >(point.x1),thickness_over_two);
}
//#####################################################################
// Function Closest_Non_Intersecting_Point
//#####################################################################
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLEX_1D<T>& point,const T thickness_over_two)
{
    RANGE<VECTOR<T,1> > thickened_box=RANGE<VECTOR<T,1> >::Bounding_Box(point.x1-2*thickness_over_two,point.x1+2*thickness_over_two);
    return Intersects(ray,RANGE<VECTOR<T,1> >(point.x1),thickness_over_two);
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,1> >&,const POINT_SIMPLEX_1D<float>&,const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,1> >&,const POINT_SIMPLEX_1D<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,1> >&,const POINT_SIMPLEX_1D<double>&,const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,1> >&,const POINT_SIMPLEX_1D<double>&,const double);
#endif
};
};
