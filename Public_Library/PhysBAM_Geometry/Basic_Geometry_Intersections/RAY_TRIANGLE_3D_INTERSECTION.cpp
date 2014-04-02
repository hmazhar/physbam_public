//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const TRIANGLE_3D<T>& triangle, const T thickness_over_two)
{
    RAY<VECTOR<T,3> > ray_temp;ray.Save_Intersection_Information(ray_temp);
    T thickness=2*thickness_over_two;

    // first check the plane of the triangle's face
    if(!INTERSECTION::Intersects(ray,static_cast<PLANE<T> >(triangle),thickness_over_two)) return false; // otherwise intersects the plane
    VECTOR<T,3> plane_point(ray.Point(ray.t_max));
    PLANE<T> edge_plane_12(VECTOR<T,3>::Cross_Product(triangle.x2-triangle.x1,triangle.normal).Normalized(),triangle.x1),edge_plane_23(VECTOR<T,3>::Cross_Product(triangle.x3-triangle.x2,triangle.normal).Normalized(),triangle.x2),
                     edge_plane_31(VECTOR<T,3>::Cross_Product(triangle.x1-triangle.x3,triangle.normal).Normalized(),triangle.x3);
    if(!edge_plane_12.Outside(plane_point,thickness) && !edge_plane_23.Outside(plane_point,thickness) && !edge_plane_31.Outside(plane_point,thickness))return true; // intersects face of triangle wedge 
    else ray.Restore_Intersection_Information(ray_temp);

    // check for intersection with the sides of the wedge
    if(edge_plane_12.Outside(plane_point,thickness_over_two) && INTERSECTION::Intersects(ray,edge_plane_12,thickness_over_two)){
        VECTOR<T,3> edge_point(ray.Point(ray.t_max));
        if(triangle.PLANE<T>::Boundary(edge_point,thickness) && !edge_plane_23.Outside(edge_point,thickness) && !edge_plane_31.Outside(edge_point,thickness)){
            ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
        else ray.Restore_Intersection_Information(ray_temp);}
    if(edge_plane_23.Outside(plane_point,thickness_over_two) && INTERSECTION::Intersects(ray,edge_plane_23,thickness_over_two)){
        VECTOR<T,3> edge_point(ray.Point(ray.t_max));
        if(triangle.PLANE<T>::Boundary(edge_point,thickness) && !edge_plane_12.Outside(edge_point,thickness) && !edge_plane_31.Outside(edge_point,thickness)){
            ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
        else ray.Restore_Intersection_Information(ray_temp);}
    if(edge_plane_31.Outside(plane_point,thickness_over_two) && INTERSECTION::Intersects(ray,edge_plane_31,thickness_over_two)){
        VECTOR<T,3> edge_point(ray.Point(ray.t_max));
        if(triangle.PLANE<T>::Boundary(edge_point,thickness) && !edge_plane_12.Outside(edge_point,thickness) && !edge_plane_23.Outside(edge_point,thickness)){
            ray.intersection_location=RAY<VECTOR<T,3> >::INTERIOR_POINT;return true;}
        else ray.Restore_Intersection_Information(ray_temp);}

    return false;
}
//#####################################################################
// Function Lazy_Intersects
//#####################################################################
template<class T> bool Lazy_Intersects(RAY<VECTOR<T,3> >& ray,const TRIANGLE_3D<T>& triangle, const T thickness_over_two)
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate_id=ray.aggregate_id;
    if(!INTERSECTION::Intersects(ray,static_cast<PLANE<T> >(triangle))) return false; // otherwise intersects the plane
    if(triangle.Lazy_Planar_Point_Inside_Triangle(ray.Point(ray.t_max))) return true; // intersects the face of the triangle 
    else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate_id;return false;} // reset ray
}
//#####################################################################
// Function Closest_Non_Intersecting_Point
//#####################################################################
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,3> >& ray,const TRIANGLE_3D<T>& triangle, const T thickness_over_two)
{
    RAY<VECTOR<T,3> > ray_temp;ray.Save_Intersection_Information(ray_temp);
    if(!INTERSECTION::Intersects(ray,triangle,thickness_over_two)) return false;
    else if(ray.intersection_location==RAY<VECTOR<T,3> >::START_POINT) return true;
    else ray.Restore_Intersection_Information(ray_temp);

    // TODO: Save having to re-generate all the planes...
    T thickness=2*thickness_over_two;
    VECTOR<T,3> normal_times_thickness=triangle.normal*thickness;
    TRIANGLE_3D<T> top_triangle(triangle.x1+normal_times_thickness,triangle.x2+normal_times_thickness,triangle.x3+normal_times_thickness);
    TRIANGLE_3D<T> bottom_triangle(triangle.x1-normal_times_thickness,triangle.x2-normal_times_thickness,triangle.x3-normal_times_thickness);
    PLANE<T> edge_plane_12(VECTOR<T,3>::Cross_Product(triangle.x2-triangle.x1,triangle.normal).Normalized(),triangle.x1);edge_plane_12.x1+=edge_plane_12.normal*thickness;
    PLANE<T> edge_plane_23(VECTOR<T,3>::Cross_Product(triangle.x3-triangle.x2,triangle.normal).Normalized(),triangle.x2);edge_plane_23.x1+=edge_plane_23.normal*thickness;
    PLANE<T> edge_plane_31(VECTOR<T,3>::Cross_Product(triangle.x1-triangle.x3,triangle.normal).Normalized(),triangle.x3);edge_plane_31.x1+=edge_plane_31.normal*thickness;
    bool found_intersection=false;
    if(INTERSECTION::Intersects(ray,top_triangle,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Intersects(ray,bottom_triangle,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Rectangle_Intersects(ray,edge_plane_12,top_triangle,bottom_triangle,edge_plane_23,edge_plane_31,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Rectangle_Intersects(ray,edge_plane_23,top_triangle,bottom_triangle,edge_plane_12,edge_plane_31,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Rectangle_Intersects(ray,edge_plane_31,top_triangle,bottom_triangle,edge_plane_12,edge_plane_23,thickness_over_two)) found_intersection=true;
    return found_intersection;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const TRIANGLE_3D<float>&,const float);
template bool Lazy_Intersects(RAY<VECTOR<float,3> >&,const TRIANGLE_3D<float>&,const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,3> >&,const TRIANGLE_3D<float>&, const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(RAY<VECTOR<double,3> >&,const TRIANGLE_3D<double>&,const double);
template bool Lazy_Intersects(RAY<VECTOR<double,3> >&,const TRIANGLE_3D<double>&,const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,3> >&,const TRIANGLE_3D<double>&, const double);
#endif
};
};
