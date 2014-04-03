//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __BOX_POLYGON_INTERSECTION_AREA__
#define __BOX_POLYGON_INTERSECTION_AREA__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class TV> typename TV::SCALAR Intersection_Area(const RANGE<TV>& box, const POLYGON<TV>& polygon);
//#####################################################################
};
};
#endif
