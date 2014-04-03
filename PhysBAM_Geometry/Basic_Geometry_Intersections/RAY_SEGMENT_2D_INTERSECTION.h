//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_SEGMENT_2D_INTERSECTION__
#define __RAY_SEGMENT_2D_INTERSECTION__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,2> >& ray,const SEGMENT_2D<T>& segment,const T thickness_over_two=0);
template<class T> bool Fuzzy_Intersects(RAY<VECTOR<T,2> >& ray,const SEGMENT_2D<T>& segment,const T thickness_over_two=0);
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,2> >& ray,const SEGMENT_2D<T>& segment,const T thickness_over_two=0);

template<class T> bool Intersection_X_Segment(RAY<VECTOR<T,2> >& ray,const T x1,const T x2,const T y,const T thickness_over_two);
template<class T> bool Intersection_Y_Segment(RAY<VECTOR<T,2> >& ray,const T x,const T y1,const T y2,const T thickness_over_two);
//#####################################################################
};
};
#endif
