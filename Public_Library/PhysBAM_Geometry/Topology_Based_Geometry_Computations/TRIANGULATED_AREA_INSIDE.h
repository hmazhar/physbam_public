//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_AREA_INSIDE__
#define __TRIANGULATED_AREA_INSIDE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T>
int Inside(const TRIANGULATED_AREA<T>& ta,const VECTOR<T,2>& location,const T thickness_over_two)
{
    PHYSBAM_ASSERT(ta.bounding_box && ta.hierarchy);
    if(ta.bounding_box->Outside(location,thickness_over_two)) return 0;
    if(ta.hierarchy->box_hierarchy(ta.hierarchy->root).Outside(location,thickness_over_two)) return 0;
    ARRAY<int> triangles_to_check;ta.hierarchy->Intersection_List(location,triangles_to_check,thickness_over_two);
    for(int l=1;l<=triangles_to_check.m;l++){
        int t=triangles_to_check(l);int i,j,k;ta.mesh.elements(t).Get(i,j,k);
        if(!TRIANGLE_2D<T>::Outside(location,ta.particles.X(i),ta.particles.X(j),ta.particles.X(k),thickness_over_two)) return t;}
    return 0;
}
template<class T>
bool Inside_Any_Simplex(const TRIANGULATED_AREA<T>& ta,const VECTOR<T,2>& location,int& triangle_id,const T thickness_over_two)
{
    assert(ta.hierarchy);
    ARRAY<int> nearby_triangles;ta.hierarchy->Intersection_List(location,nearby_triangles,thickness_over_two);
    for(int k=1;k<=nearby_triangles.m;k++){
        const TRIANGLE_2D<T>& triangle=ta.Get_Element(nearby_triangles(k));
        if(triangle.Inside(location,thickness_over_two)){triangle_id=nearby_triangles(k);return true;}}
    return false;
}
}
}
#endif
