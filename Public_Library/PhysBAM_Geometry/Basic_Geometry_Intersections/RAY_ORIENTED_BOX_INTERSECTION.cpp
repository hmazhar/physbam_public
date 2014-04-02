//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_ORIENTED_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Fuzzy_Intersects(RAY<VECTOR<T,2> >& ray,const ORIENTED_BOX<VECTOR<T,2> >& box,const T segment_intersect_epsilon)
{
    bool intersects=false;
    for(int i=1;i<=2;++i) if(INTERSECTION::Fuzzy_Intersects(ray,SEGMENT_2D<T>(box.corner,box.corner+box.edges.Column(i)),segment_intersect_epsilon)) intersects=true;
    for(int i=1;i<=2;++i) if(INTERSECTION::Fuzzy_Intersects(ray,SEGMENT_2D<T>(box.corner+box.edges.Column(i),box.corner+box.edges.Column_Sum()),segment_intersect_epsilon)) intersects=true;
    return intersects;
}
//#####################################################################
template bool Fuzzy_Intersects(RAY<VECTOR<float,2> >&,const ORIENTED_BOX<VECTOR<float,2> >&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Fuzzy_Intersects(RAY<VECTOR<double,2> >&,const ORIENTED_BOX<VECTOR<double,2> >&,const double);
#endif
};
};
