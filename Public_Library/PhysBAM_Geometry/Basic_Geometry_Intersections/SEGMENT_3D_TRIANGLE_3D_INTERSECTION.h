//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __SEGMENT_3D_TRIANGLE_3D_INTERSECTION__
#define __SEGMENT_3D_TRIANGLE_3D_INTERSECTION__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(const SEGMENT_3D<T>& segment,const TRIANGLE_3D<T>& triangle,const T thickness_over_two=0);
template<class T> bool Intersects(const SEGMENT_3D<T>& segment,const TRIANGLE_3D<T>& triangle,T& a,VECTOR<T,3>& weights,const T thickness_over_two=0);
//#####################################################################
};
};
#endif
