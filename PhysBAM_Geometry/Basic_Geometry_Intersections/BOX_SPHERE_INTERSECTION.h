//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __BOX_SPHERE_INTERSECTION__
#define __BOX_SPHERE_INTERSECTION__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T,class TV> bool Intersects(const RANGE<TV>&,const SPHERE<TV>&,const T& thickness_over_two=0);
//#####################################################################
};
};
#endif
