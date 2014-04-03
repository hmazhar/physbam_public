//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_INSIDE_USING_RAY_TEST__
#define __TRIANGULATED_SURFACE_INSIDE_USING_RAY_TEST__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;
template<class TV> class RAY;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T> bool
Inside_Using_Ray_Test(const TRIANGULATED_SURFACE<T>& ts,RAY<VECTOR<T,3> >& ray,const T thickness_over_two)
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(ts.mesh.adjacent_elements && ts.triangle_list);
    bool inside=false;
    if(INTERSECTION::Intersects(ray,ts,thickness_over_two) && ray.t_max > 0){ // otherwise missed the object or on the boundary
        TV point=ray.Point(ray.t_max); // point is inside if and only if location is inside
        T thickness=2*thickness_over_two;
        TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(ray.aggregate_id);
        int region_id,region=triangle.Region(point,region_id,thickness);
        if(region == 1){ // vertex
            if(ts.Signed_Solid_Angle_Of_Triangle_Web(point,ts.mesh.elements(ray.aggregate_id)(region_id)) > 0) inside=true;} 
        else if(region == 2) { // edge
            int node1,node2,node3,neighbor=0;ts.mesh.elements(ray.aggregate_id).Get(node1,node2,node3);
            if(region_id==1) neighbor=ts.mesh.Adjacent_Triangle(ray.aggregate_id,node1,node2);
            else if(region_id==2) neighbor=ts.mesh.Adjacent_Triangle(ray.aggregate_id,node2,node3);
            else neighbor=ts.mesh.Adjacent_Triangle(ray.aggregate_id,node3,node1);
            if(neighbor==0){if(triangle.Lazy_Inside(ray.endpoint)) inside=true;}
            else{
                TRIANGLE_3D<T>& triangle2=(*ts.triangle_list)(neighbor);
                int convex=0;
                if(region_id == 1){if(TV::Dot_Product(triangle2.normal,triangle.x3-triangle2.x1) >= 0) convex=1;}
                else if(region_id == 2){if(TV::Dot_Product(triangle2.normal,triangle.x1-triangle2.x1) >= 0) convex=1;}
                else if(region_id == 3){if(TV::Dot_Product(triangle2.normal,triangle.x2-triangle2.x1) >= 0) convex=1;}
                if(convex){if(triangle.Lazy_Inside(ray.endpoint) && triangle2.Lazy_Inside(ray.endpoint)) inside=true;} // inside both - can use location or point
                else{if(triangle.Lazy_Inside(ray.endpoint) || triangle2.Lazy_Inside(ray.endpoint)) inside=true;}}} // inside either - can use location or point
        else{if(triangle.Lazy_Inside(ray.endpoint)) inside=true;}} // region=3 - face - can use location or point
    return inside;
}
}
}
#endif
