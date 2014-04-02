//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_PLANE_INTERSECTION__
#define __RAY_PLANE_INTERSECTION__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane, const T thickness_over_two,const T distance,const T rate_of_approach);
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane, const T thickness_over_two=0);
template<class T> bool Lazy_Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane);
template<class T> bool Rectangle_Intersects(RAY<VECTOR<T,3> >& ray,const PLANE<T>& plane1,const PLANE<T>& plane2,const PLANE<T>& plane3,const PLANE<T>& plane4,const PLANE<T>& plane5,const T thickness_over_two=0);
//#####################################################################
}
}
#endif
