//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Andrew Selle, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_3D
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Change_Size
//#####################################################################
// Enlarges the triangle by pushing out the triangle edges by distance 'delta' orthogonally to the edges.
// This keeps the incenter fixed.  If the triangle is degenerate, it will not be changed.
template<class T> void TRIANGLE_3D<T>::
Change_Size(const T delta)
{
    VECTOR<T,3> edge_lengths((x3-x2).Magnitude(),(x1-x3).Magnitude(),(x2-x1).Magnitude());
    T perimeter=edge_lengths.Sum(),area=Area();
    if(!perimeter || !area) return; // don't know which direction to enlarge a degenerate triangle, so do nothing
    T scale=1+delta*(T).5*perimeter/area;
    VECTOR<T,3> incenter=Point_From_Barycentric_Coordinates(edge_lengths/perimeter);
    x1=incenter+(x1-incenter)*scale;x2=incenter+(x2-incenter)*scale;x3=incenter+(x3-incenter)*scale;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Inside(const TV& point,const T thickness_over_two) const
{
    return Point_Inside_Triangle(point,thickness_over_two);
}
//#####################################################################
// Function Point_Inside_Triangle
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Point_Inside_Triangle(const VECTOR<T,3>& point,const T thickness_over_2) const
{
    return PLANE<T>::Boundary(point,thickness_over_2)&&Planar_Point_Inside_Triangle(point,thickness_over_2);
}
//#####################################################################
// Function Planar_Point_Inside_Triangle
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Planar_Point_Inside_Triangle(const VECTOR<T,3>& point,const T thickness_over_2) const
{
    PLANE<T> edge_plane(VECTOR<T,3>::Cross_Product(x2-x1,normal).Normalized(),x1);if(edge_plane.Outside(point,thickness_over_2)) return false;
    edge_plane.normal=VECTOR<T,3>::Cross_Product(x1-x3,normal).Normalized();if(edge_plane.Outside(point,thickness_over_2)) return false;
    edge_plane.normal=VECTOR<T,3>::Cross_Product(x3-x2,normal).Normalized();edge_plane.x1=x2;if(edge_plane.Outside(point,thickness_over_2)) return false;
    return true;
}
//#####################################################################
// Function Lazy_Planar_Point_Inside_Triangle
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Lazy_Planar_Point_Inside_Triangle(const VECTOR<T,3>& point) const
{
    VECTOR<T,3> edge_normal_1=VECTOR<T,3>::Cross_Product(x2-x1,normal),point_minus_x1=point-x1;if(VECTOR<T,3>::Dot_Product(edge_normal_1,point_minus_x1) > 0) return false;
    VECTOR<T,3> edge_normal_2=VECTOR<T,3>::Cross_Product(x1-x3,normal);if(VECTOR<T,3>::Dot_Product(edge_normal_2,point_minus_x1) > 0) return false;
    edge_normal_1+=edge_normal_2;if(VECTOR<T,3>::Dot_Product(edge_normal_1,point-x2) < 0) return false; // this equals x2-x3 (== -edge_normal_3)
    return true;
}
//#####################################################################
// Function Minimum_Edge_Length
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Minimum_Edge_Length() const
{      
    return min((x2-x1).Magnitude(),(x3-x2).Magnitude(),(x1-x3).Magnitude());
}
//#####################################################################
// Function Maximum_Edge_Length
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Maximum_Edge_Length() const
{      
    return max((x2-x1).Magnitude(),(x3-x2).Magnitude(),(x1-x3).Magnitude());
}
//#####################################################################
// Function Region
//#####################################################################
// returns the region one is near in priority of 1=vertex, 2=edge, 3=face based on distance, region_id differentiates which point or edge
template<class T> int TRIANGLE_3D<T>::
Region(const VECTOR<T,3>& location,int& region_id,const T distance) const
{
    SEGMENT_3D<T> segment(x1,x2);T d1=segment.Distance_From_Point_To_Segment(location);
    segment.x1=x3;T d2=segment.Distance_From_Point_To_Segment(location);
    segment.x2=x1;T d3=segment.Distance_From_Point_To_Segment(location);
    if(d1 > distance && d2 > distance && d3 > distance) return 3; // face is closest
    else if(d1 <= distance){
        if(d2 <= distance){region_id=2;return 1;} // vertex 2
        else if(d3 <= distance){region_id=1;return 1;} // vertex 1
        else{region_id=1;return 2;}} // edge 1
    else if(d2 <= distance){
        if(d3 <= distance){region_id=3;return 1;} // vertex 3
        else{region_id=2;return 2;}} // edge 2
    else{region_id=3;return 2;} // edge 3
}
//#####################################################################
// Function Closest_Point
//#####################################################################
template<class T> VECTOR<T,3> TRIANGLE_3D<T>::
Closest_Point(const VECTOR<T,3>& location,VECTOR<T,3>& weights) const
{
    weights=Barycentric_Coordinates(location);
    // project closest point to the triangle if it's not already inside it
    if(weights.x<0){
        T a23=SEGMENT_3D<T>::Interpolation_Fraction(location,x2,x3); // Check edge x2--x3
        if(a23<0){
            if(weights.z<0){ // Closest point is on edge x1--x2
                T a12=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,x1,x2),0,1);weights=VECTOR<T,3>(1-a12,a12,0);return weights.x*x1+weights.y*x2;}
            else{weights=VECTOR<T,3>(0,1,0);return x2;}} // Closest point is x2
        else if(a23>1){
            if(weights.y<0){ // Closest point is on edge x1--x3
                T a13=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,x1,x3),0,1);weights=VECTOR<T,3>(1-a13,0,a13);return weights.x*x1+weights.z*x3;}
            else{weights=VECTOR<T,3>(0,0,1);return x3;}} // Closest point is x3
        else{weights=VECTOR<T,3>(0,1-a23,a23);return weights.y*x2+weights.z*x3;}} // Closest point is on edge x2--x3
    else if(weights.y<0){
        T a13=SEGMENT_3D<T>::Interpolation_Fraction(location,x1,x3); // Check edge x1--x3
        if(a13<0){
            if(weights.z<0){ // Closest point is on edge x1--x2
                T a12=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,x1,x2),0,1);weights=VECTOR<T,3>(1-a12,a12,0);return weights.x*x1+weights.y*x2;}
            else{weights=VECTOR<T,3>(1,0,0);return x1;}} // Closest point is x1
        else if(a13>1){weights=VECTOR<T,3>(0,0,1);return x3;} // Closest point is x3
        else{weights=VECTOR<T,3>(1-a13,0,a13);return weights.x*x1+weights.z*x3;}} // Closest point is on edge x1--x3
    else if(weights.z<0){ // Closest point is on edge x1--x2
        T a12=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,x1,x2),0,1);weights=VECTOR<T,3>(1-a12,a12,0);return weights.x*x1+weights.y*x2;}
    return weights.x*x1+weights.y*x2+weights.z*x3; // Point is interior to the triangle
}
//#####################################################################
// Function Distance_To_Triangle
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Distance_To_Triangle(const VECTOR<T,3>& location) const
{   
    VECTOR<T,3> weights,projected_point;
    projected_point=Closest_Point(location,weights);return (location-projected_point).Magnitude();
}
//#####################################################################
// Function Minimum_Angle
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Minimum_Angle() const
{
    VECTOR<T,3> s1=(x1-x2).Normalized(),s2=(x2-x3).Normalized(),s3=(x3-x1).Normalized();
    return acos(max(VECTOR<T,3>::Dot_Product(s1,-s2),VECTOR<T,3>::Dot_Product(-s1,s3),VECTOR<T,3>::Dot_Product(s2,-s3)));
}
//#####################################################################
// Function Maximum_Angle
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Maximum_Angle() const
{
    VECTOR<T,3> s1=(x1-x2).Normalized(),s2=(x2-x3).Normalized(),s3=(x3-x1).Normalized();
    return acos(min(VECTOR<T,3>::Dot_Product(s1,-s2),VECTOR<T,3>::Dot_Product(-s1,s3),VECTOR<T,3>::Dot_Product(s2,-s3)));
}
//#####################################################################
// Function Signed_Solid_Angle
//#####################################################################
// positive for normals that point away from the center - not reliable if center is too close to the triangle face
template<class T> T TRIANGLE_3D<T>::
Signed_Solid_Angle(const VECTOR<T,3>& center) const
{
    VECTOR<T,3> r=(x1-center).Normalized(),u=x2-x1,v=x3-x1;u-=VECTOR<T,3>::Dot_Product(u,r)*r;v-=VECTOR<T,3>::Dot_Product(v,r)*r;
    T solid_angle=-(T)pi+VECTOR<T,3>::Angle_Between(u,v);
    r=(x2-center).Normalized();u=x1-x2,v=x3-x2;u-=VECTOR<T,3>::Dot_Product(u,r)*r;v-=VECTOR<T,3>::Dot_Product(v,r)*r;
    solid_angle+=VECTOR<T,3>::Angle_Between(u,v);
    r=(x3-center).Normalized();u=x1-x3,v=x2-x3;u-=VECTOR<T,3>::Dot_Product(u,r)*r;v-=VECTOR<T,3>::Dot_Product(v,r)*r;
    solid_angle+=VECTOR<T,3>::Angle_Between(u,v);
    solid_angle=max(T(0),min((T)(2*pi),solid_angle));
    if(VECTOR<T,3>::Dot_Product(r,normal) < 0) solid_angle*=(-1);
    return solid_angle;
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
// outputs unsigned distance
template<class T> bool TRIANGLE_3D<T>::
Point_Face_Interaction(const VECTOR<T,3>& x,const T interaction_distance,const bool allow_negative_weights,T& distance) const                       
{      
    distance=VECTOR<T,3>::Dot_Product(x-x1,normal);
    return abs(distance)<=interaction_distance && Planar_Point_Inside_Triangle(x,allow_negative_weights?interaction_distance:0);
}
//#####################################################################
// Function Point_Face_Interaction_Data
//#####################################################################
template<class T> void TRIANGLE_3D<T>::
Point_Face_Interaction_Data(const VECTOR<T,3>& x,T& distance,VECTOR<T,3>& interaction_normal,VECTOR<T,3>& weights,const bool perform_attractions) const                       
{
    interaction_normal=normal;weights=Barycentric_Coordinates(x);
    if(!perform_attractions && distance<0){distance*=-1;interaction_normal*=-1;} // distance > 0, interaction_normal points from the triangle to the point
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Point_Face_Interaction(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T interaction_distance,T& distance,
    VECTOR<T,3>& interaction_normal,VECTOR<T,3>& weights,T& relative_speed,const bool allow_negative_weights,const bool exit_early) const
{      
    if(!Point_Face_Interaction(x,interaction_distance,allow_negative_weights,distance)) return false;
    if(!exit_early){
        Point_Face_Interaction_Data(x,distance,interaction_normal,weights,false);
        relative_speed=TV::Dot_Product(v-(weights.x*v1+weights.y*v2+weights.z*v3),interaction_normal);} // relative speed is in the normal direction
    return true;
}
//#####################################################################
// Function Robust_Point_Triangle_Collision
//#####################################################################
template<class T> POINT_SIMPLEX_COLLISION_TYPE TRIANGLE_3D<T>::
Robust_Point_Triangle_Collision(const TRIANGLE_3D<T>& initial_triangle,const TRIANGLE_3D<T>& final_triangle,const VECTOR<T,3>& x,const VECTOR<T,3>& final_x,const T dt,
    const T collision_thickness,T& collision_time,VECTOR<T,3>& normal,VECTOR<T,3>& weights,T& relative_speed)
{
    return CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Triangle_Collision(initial_triangle,final_triangle,x,final_x,dt,collision_thickness,collision_time,
        normal,weights,relative_speed);
}
//#####################################################################
// Function Point_Triangle_Collision
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Point_Face_Collision(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T dt,const T collision_thickness,
    T& collision_time,TV& normal,TV& weights,T& relative_speed,const bool exit_early) const
{
    return CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Point_Face_Collision(*this,x,v,v1,v2,v3,dt,collision_thickness,collision_time,normal,weights,relative_speed,exit_early);
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void TRIANGLE_3D<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_3D<T> >& clipped_simplices) const
{
    // cut with all sides of box
    clipped_simplices.Remove_All();
    clipped_simplices.Append(*this);
    for(int axis=1;axis<=TV::dimension;axis++){
        for(int i=clipped_simplices.m;i>=1;i--){
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),PLANE<T>(-TV::Axis_Vector(axis),box.min_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}
        for(int i=clipped_simplices.m;i>=1;i--){
            // TODO: make this more efficient by not removing a triangle that is fully inside
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),PLANE<T>(TV::Axis_Vector(axis),box.max_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}}
}
//#####################################################################
// Function Cut_With_Hyperplane_And_Discard_Outside_Simplices
//#####################################################################
template<class T> void TRIANGLE_3D<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_3D<T>& triangle,const PLANE<T>& cutting_plane,ARRAY<TRIANGLE_3D<T> >& negative_triangles)
{
    VECTOR<T,3> phi_nodes;
    VECTOR<VECTOR<T,3>,3> X_nodes;X_nodes(1)=triangle.x1;X_nodes(2)=triangle.x2;X_nodes(3)=triangle.x3;
    for(int i=1;i<=3;i++){phi_nodes(i)=cutting_plane.Signed_Distance(X_nodes(i));}

    // left simplices are in the negative halfspace, right simplices in the positive halfspace
    int positive_count=0,single_node_sign;
    for(int i=1;i<=3;i++) if(phi_nodes(i)>0) positive_count++;
    switch(positive_count){
        case 0: // in negative halfspace
            negative_triangles.Append(triangle);break;
        case 1:
        case 2:
            single_node_sign=positive_count==1?1:-1;
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            for(int i=1;i<=3;i++)if(LEVELSET_UTILITIES<T>::Sign(phi_nodes(i))==single_node_sign){
                VECTOR<VECTOR<T,3>,2> interface_locations;int index=i%3+1;
                VECTOR<int,2> other_locations;
                for(int j=1;j<=2;j++,index=index%3+1){
                    other_locations(j)=index;
                    interface_locations(j)=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes(i),X_nodes(index),LEVELSET_UTILITIES<T>::Theta(phi_nodes(i),phi_nodes(index)));}
                if(positive_count==1){ // add two triangles to negative triangles
                    negative_triangles.Append(TRIANGLE_3D<T>(interface_locations(1),X_nodes(other_locations(1)),X_nodes(other_locations(2))));
                    negative_triangles.Append(TRIANGLE_3D<T>(interface_locations(1),X_nodes(other_locations(2)),interface_locations(2)));}
                else // add triangle to negative_triangles
                    negative_triangles.Append(TRIANGLE_3D<T>(X_nodes(i),interface_locations(1),interface_locations(2)));
                return;}
        case 3: // in positive halfspace
            break;}
}
//#####################################################################
template class TRIANGLE_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_3D<double>;
#endif
