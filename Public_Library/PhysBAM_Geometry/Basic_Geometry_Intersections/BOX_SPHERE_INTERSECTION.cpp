//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_SPHERE_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
// From http://www.acm.org/tog/GraphicsGems/gems/BoxSphere.c
template<class T,class TV> bool Intersects(const RANGE<TV>& box,const SPHERE<TV>& sphere,const T& thickness_over_two)
{
    T dmin=0;TV min_diff=sphere.center-box.min_corner+thickness_over_two,max_diff=sphere.center-box.max_corner-thickness_over_two;
    for(int i=1;i<=TV::dimension;i++){if(min_diff(i)<0) dmin+=sqr(min_diff(i));else if(max_diff(i)>0) dmin+=sqr(max_diff(i));}
    return dmin<=sqr(sphere.radius);
}
//#####################################################################
template bool Intersects(const RANGE<VECTOR<float,1> >&,const SPHERE<VECTOR<float,1> >&,const float&);
template bool Intersects(const RANGE<VECTOR<float,2> >&,const SPHERE<VECTOR<float,2> >&,const float&);
template bool Intersects(const RANGE<VECTOR<float,3> >&,const SPHERE<VECTOR<float,3> >&,const float&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const RANGE<VECTOR<double,1> >&,const SPHERE<VECTOR<double,1> >&,const double&);
template bool Intersects(const RANGE<VECTOR<double,2> >&,const SPHERE<VECTOR<double,2> >&,const double&);
template bool Intersects(const RANGE<VECTOR<double,3> >&,const SPHERE<VECTOR<double,3> >&,const double&);
#endif
};
};
