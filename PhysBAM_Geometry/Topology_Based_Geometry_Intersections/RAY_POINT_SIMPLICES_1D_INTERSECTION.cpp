//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_POINT_SIMPLEX_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_POINT_SIMPLICES_1D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLICES_1D<T>& simplices, const T thickness_over_two)
{
    bool hit=false;
    for(int i=1;i<=simplices.mesh.elements.m;++i){
        POINT_SIMPLEX_1D<T> point_simplex=simplices.point_simplex_list?(*simplices.point_simplex_list)(i):simplices.Get_Element(i);
        if(INTERSECTION::Intersects(ray,point_simplex,thickness_over_two)){ray.aggregate_id=i;hit=true;}}
    return hit;
}
//#####################################################################
// Function Closest_Non_Intersecting_Point
//#####################################################################
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLICES_1D<T>& simplices,const T thickness_over_two)
{
    bool hit=false;
    for(int i=1;i<=simplices.mesh.elements.m;++i){
        POINT_SIMPLEX_1D<T> point_simplex=simplices.point_simplex_list?(*simplices.point_simplex_list)(i):simplices.Get_Element(i);
        if(INTERSECTION::Closest_Non_Intersecting_Point(ray,point_simplex,thickness_over_two)){ray.aggregate_id=i;hit=true;}}
    return hit;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,1> >&,const POINT_SIMPLICES_1D<float>&, const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,1> >&,const POINT_SIMPLICES_1D<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,1> >&,const POINT_SIMPLICES_1D<double>&, const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,1> >&,const POINT_SIMPLICES_1D<double>&,const double);
#endif
};
};
