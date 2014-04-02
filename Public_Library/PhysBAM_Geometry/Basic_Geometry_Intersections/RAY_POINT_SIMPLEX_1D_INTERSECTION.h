//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_POINT_SIMPLEX_1D_INTERSECTION__
#define __RAY_POINT_SIMPLEX_1D_INTERSECTION__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLEX_1D<T>& plane,const T thickness_over_two=0);
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLEX_1D<T>& segment,const T thickness_over_two=0);
//#####################################################################
};
};
#endif
