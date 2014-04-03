//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_SURFACE__
#define __TRIANGULATED_SURFACE_SURFACE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>

namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class TV> class RAY;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T>
VECTOR<T,3> Surface(const TRIANGULATED_SURFACE<T>& ts,const VECTOR<T,3>& location,const T max_depth,const T thickness_over_2,int* closest_triangle,T* distance)
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(ts.triangle_list);
    TV point;
    
    if(max_depth){
        PHYSBAM_ASSERT(ts.hierarchy);
        RANGE<TV> box(location.x-max_depth,location.x+max_depth,location.y-max_depth,location.y+max_depth,location.z-max_depth,location.z+max_depth);
        ARRAY<int> nearby_triangles;ts.hierarchy->Intersection_List(box,nearby_triangles);
        if(!nearby_triangles.m){ // grab any point assuming far from the interface
            RAY<TV> ray(location,ts.particles.X(ts.mesh.elements(1)(1))-location);if(closest_triangle) *closest_triangle=1;
            if(INTERSECTION::Intersects(ray,ts,thickness_over_2)){if(distance) *distance=ray.t_max;return ray.Point(ray.t_max);}
            else{if(distance) *distance=(ts.particles.X(ts.mesh.elements(1)(1))-location).Magnitude();return ts.particles.X(ts.mesh.elements(1)(1));}} 
        else{
            TV weights;
            {TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(nearby_triangles(1));point=triangle.Closest_Point(location,weights);}
            T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=nearby_triangles(1);
            for(int k=2;k<=nearby_triangles.m;k++){
                TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(nearby_triangles(k));
                TV new_point=triangle.Closest_Point(location,weights);
                T new_distance=(location-new_point).Magnitude_Squared();
                if(new_distance < distance_temp){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=nearby_triangles(k);}}
            if(distance) *distance=sqrt(distance_temp);return point;}}

    // slow method
    TV weights;
    {TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(1);point=triangle.Closest_Point(location,weights);}
    T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=1;
    for(int k=2;k<=ts.mesh.elements.m;k++){
        TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(k);
        TV new_point=triangle.Closest_Point(location,weights);
        T new_distance=(location-new_point).Magnitude_Squared();
        if(new_distance < distance_temp){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=k;}}
    if(distance) *distance=sqrt(distance_temp);return point;
}
template<class T>
VECTOR<T,3> Oriented_Surface(const TRIANGULATED_SURFACE<T>& ts,const VECTOR<T,3>& location,const VECTOR<T,3>& normal,const T max_depth,const T thickness_over_2,int* closest_triangle,T* distance)
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(ts.triangle_list);
    TV point;
    
    if(max_depth){
        PHYSBAM_ASSERT(ts.hierarchy);
        RANGE<TV> box(location.x-max_depth,location.x+max_depth,location.y-max_depth,location.y+max_depth,location.z-max_depth,location.z+max_depth);
        ARRAY<int> nearby_triangles;ts.hierarchy->Intersection_List(box,nearby_triangles);
        if(!nearby_triangles.m){ // grab any point ignoring normal assuming far from the interface
            RAY<TV> ray(location,ts.particles.X(ts.mesh.elements(1)(1))-location);if(closest_triangle) *closest_triangle=1;
            if(INTERSECTION::Intersects(ray,ts,thickness_over_2)){if(distance) *distance=ray.t_max;return ray.Point(ray.t_max);}
            else{if(distance) *distance=(ts.particles.X(ts.mesh.elements(1)(1))-location).Magnitude();return ts.particles.X(ts.mesh.elements(1)(1));}} 
        else{
            TV weights;
            {TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(nearby_triangles(1));point=triangle.Closest_Point(location,weights);}
            T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=nearby_triangles(1);
            for(int k=2;k<=nearby_triangles.m;k++){
                TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(nearby_triangles(k));
                TV new_point=triangle.Closest_Point(location,weights);
                T new_distance=(location-new_point).Magnitude_Squared();
                if(new_distance<distance_temp&&TV::Dot_Product(normal,triangle.normal)>0){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=nearby_triangles(k);}}
            if(distance) *distance=sqrt(distance_temp);return point;}}

    // slow method
    TV weights;
    {TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(1);point=triangle.Closest_Point(location,weights);}
    T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=1;
    for(int k=2;k<=ts.mesh.elements.m;k++){
        TRIANGLE_3D<T>& triangle=(*ts.triangle_list)(k);
        TV new_point=triangle.Closest_Point(location,weights);
        T new_distance=(location-new_point).Magnitude_Squared();
        if(new_distance<distance_temp&&TV::Dot_Product(normal,triangle.normal)>0){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=k;}}
    if(distance) *distance=sqrt(distance_temp);return point;
}

}
}
#endif
