//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const RANGE<VECTOR<T,3> >& box,const TRIANGLE_3D<T>& triangle,const T thickness_over_two)
{
    typedef VECTOR<T,3> TV;
    TV c=box.Center(),h=(T).5*box.Edge_Lengths()+thickness_over_two,x1=triangle.x1-c,x2=triangle.x2-c,x3=triangle.x3-c;
    VECTOR<TV,3> cp(TV::Cross_Product(x1,x2),TV::Cross_Product(x2,x3),TV::Cross_Product(x3,x1));
    TV sum(cp.Sum()),minsum=TV::Componentwise_Min(sum,TV()),maxsum=TV::Componentwise_Max(sum,TV());
    VECTOR<TV,3> diff=abs(VECTOR<TV,3>(x1-x2,x2-x3,x3-x1));

    for(int i=1;i<=3;i++) for(int j=1;j<=3;j++){int k=j+1-3*(j==3),m=6-k-j;
        T r=h(k)*diff(i)(m)+h(m)*diff(i)(k);if(minsum(j)-cp(i)(j)>r || maxsum(j)-cp(i)(j)<-r) return false;}

    if(!box.Intersection(triangle.Bounding_Box(),thickness_over_two)) return false;
    if(!Intersects(box,static_cast<const PLANE<T>&>(triangle),thickness_over_two)) return false;

    return true;
}
//#####################################################################
template bool Intersects(const RANGE<VECTOR<float,3> >&,const TRIANGLE_3D<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const RANGE<VECTOR<double,3> >&,const TRIANGLE_3D<double>&,const double);
#endif
};
};
