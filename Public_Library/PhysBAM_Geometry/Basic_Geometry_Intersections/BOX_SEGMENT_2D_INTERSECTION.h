//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __BOX_SEGMENT_2D_INTERSECTION__
#define __BOX_SEGMENT_2D_INTERSECTION__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(const RANGE<VECTOR<T,2> >& box,const SEGMENT_2D<T>& plane,const T thickness_over_two=0);
//#####################################################################
};
};
#endif
