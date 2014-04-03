//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_INSIDE__
#define __TRIANGULATED_SURFACE_INSIDE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>

namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class TV> class RAY;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T>
bool Inside(const TRIANGULATED_SURFACE<T>& ts,const VECTOR<T,3>& location,const T thickness_over_two) 
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(ts.bounding_box && ts.hierarchy);
    if(ts.bounding_box->Outside(location,thickness_over_two)) return false;
    if(ts.hierarchy->box_hierarchy(ts.hierarchy->root).Outside(location,thickness_over_two)) return false;
    RAY<TV> ray(location,TV(0,0,1),true);
    return ts.Inside_Using_Ray_Test(ray,thickness_over_two);
}
template<class T>
bool Outside(const TRIANGULATED_SURFACE<T>& ts,const VECTOR<T,3>& location,const T thickness_over_two)
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(ts.bounding_box && ts.hierarchy && ts.mesh.adjacent_elements && ts.triangle_list);
    if(ts.bounding_box->Outside(location,thickness_over_two)) return true;
    if(ts.hierarchy->box_hierarchy(ts.hierarchy->root).Outside(location,thickness_over_two)) return true;
       
    bool outside=false;
    RAY<TV> ray(location,TV(0,0,1));
    if(!INTERSECTION::Intersects(ray,ts,thickness_over_two)) outside=true; // missed the object, outside
    else if(ray.t_max > 0){ // not in boundary region
        TV point=ray.Point(ray.t_max); // point is inside if and only if location is inside
        T thickness=2*thickness_over_two;
        TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(ray.aggregate_id);
        int region_id,region=triangle.Region(point,region_id,thickness);
        if(region == 1){ // vertex
            if(ts.Signed_Solid_Angle_Of_Triangle_Web(point,ts.mesh.elements(ray.aggregate_id)(region_id)) < 0) outside=true;} 
        else if(region == 2){ // edge
            int node1,node2,node3,neighbor=0;ts.mesh.elements(ray.aggregate_id).Get(node1,node2,node3);
            if(region_id==1) neighbor=ts.mesh.Adjacent_Triangle(ray.aggregate_id,node1,node2);
            else if(region_id==2) neighbor=ts.mesh.Adjacent_Triangle(ray.aggregate_id,node2,node3);
            else neighbor=ts.mesh.Adjacent_Triangle(ray.aggregate_id,node3,node1);
            if(neighbor==0){if(triangle.Lazy_Outside(location)) outside=true;}
            else{
                TRIANGLE_3D<T>& triangle2=(*ts.triangle_list)(neighbor);
                int convex=0;
                if(region_id == 1){if(TV::Dot_Product(triangle2.normal,triangle.x3-triangle2.x1) >= 0) convex=1;}
                else if(region_id == 2){if(TV::Dot_Product(triangle2.normal,triangle.x1-triangle2.x1) >= 0) convex=1;}
                else if(region_id == 3){if(TV::Dot_Product(triangle2.normal,triangle.x2-triangle2.x1) >= 0) convex=1;}
                if(convex){if(triangle.Lazy_Outside(location) || triangle2.Lazy_Outside(location)) outside=true;} // outside either - can use location or point
                else{if(triangle.Lazy_Outside(location) && triangle2.Lazy_Outside(location)) outside=true;}}} // outside both - can use location or point
        else{if(triangle.Lazy_Outside(location)) outside=true;}} // region=3 - face - can use location or point
    return outside;
}
template<class T>
bool Inside_Any_Triangle(const TRIANGULATED_SURFACE<T>& ts,const VECTOR<T,3>& location,int& triangle_id,const T thickness_over_two)
{
    assert(ts.hierarchy);assert(ts.triangle_list);
    ARRAY<int> nearby_triangles;ts.hierarchy->Intersection_List(location,nearby_triangles,thickness_over_two);
    for(int k=1;k<=nearby_triangles.m;k++){
        TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(nearby_triangles(k));
        if(triangle.Point_Inside_Triangle(location,thickness_over_two)){triangle_id=nearby_triangles(k);return true;}}
    return false;
}
}
}
#endif
