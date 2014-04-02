//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_SPHERE_INTERSECTION__
#define __RAY_SPHERE_INTERSECTION__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const SPHERE<VECTOR<T,1> >& box, const T thickness=0);
template<class T> bool Intersects(RAY<VECTOR<T,2> >& ray,const SPHERE<VECTOR<T,2> >& box, const T thickness=0);
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const SPHERE<VECTOR<T,3> >& box, const T thickness=0);
//#####################################################################
}
}
#endif
