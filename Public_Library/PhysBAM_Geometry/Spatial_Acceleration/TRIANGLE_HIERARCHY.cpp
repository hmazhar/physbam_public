//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGLE_HIERARCHY<T>::
TRIANGLE_HIERARCHY(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >& particles_input,const bool update_boxes,const int triangles_per_group_input)
    :triangle_mesh(triangle_mesh_input),particles(particles_input),triangle_list(0),triangles_per_group(triangles_per_group_input)
{
    if(triangle_mesh.elements.m){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGLE_HIERARCHY<T>::
TRIANGLE_HIERARCHY(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >& particles_input,ARRAY<TRIANGLE_3D<T> >& triangle_list_input,const bool update_boxes,const int triangles_per_group_input)
    :triangle_mesh(triangle_mesh_input),particles(particles_input),triangle_list(&triangle_list_input),triangles_per_group(triangles_per_group_input)
{
    if(triangle_mesh.elements.m){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TRIANGLE_HIERARCHY<T>::
~TRIANGLE_HIERARCHY()
{}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class T> void TRIANGLE_HIERARCHY<T>::
Initialize_Hierarchy_Using_KD_Tree()
{
    KD_TREE<VECTOR<T,3> > kd_tree(false);
    ARRAY<VECTOR<T,3> > centroids(triangle_mesh.elements.m);
    for(int t=1;t<=triangle_mesh.elements.m;t++){int i,j,k;triangle_mesh.elements(t).Get(i,j,k);centroids(t)=TRIANGLE_3D<T>::Center(particles.X(i),particles.X(j),particles.X(k));}
    if(triangles_per_group){
        triangles_in_group.Clean_Memory();
        kd_tree.Create_Left_Balanced_KD_Tree_With_Grouping(centroids,triangles_in_group,triangles_per_group);
        leaves=triangles_in_group.m;parents.Resize(leaves);}
    else{
        kd_tree.Create_Left_Balanced_KD_Tree(centroids);
        leaves=triangle_mesh.elements.m;parents.Resize(leaves);}
    children.Remove_All();root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);
    assert(root==2*leaves-1);
    box_hierarchy.Resize(root);box_radius.Resize(root);
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class T> void TRIANGLE_HIERARCHY<T>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes)
{
    Calculate_Bounding_Boxes(bounding_boxes,particles.X);
}
//#####################################################################
// Function Calculate_Bounding_Boxes_Helper
//#####################################################################
template<class T> template<class T_ARRAY_TV> void TRIANGLE_HIERARCHY<T>::
Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV X) 
{
    STATIC_ASSERT((IS_SAME<TV,typename T_ARRAY_TV::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY_TV>::value));
    if(triangles_per_group) for(int k=1;k<=leaves;k++){
        if(triangles_in_group(k).m) bounding_boxes(k)=RANGE<TV>::Bounding_Box(X.Subset(triangle_mesh.elements(triangles_in_group(k)(1))));
        for(int i=2;i<=triangles_in_group(k).m;i++) bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(X.Subset(triangle_mesh.elements(triangles_in_group(k)(i))));}
    else for(int k=1;k<=leaves;k++)
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(X.Subset(triangle_mesh.elements(k)));
}
//#####################################################################
// Function Calculate_Bounding_Boxes_Helper
//#####################################################################
template<class T> template<class T_ARRAY_TV> void TRIANGLE_HIERARCHY<T>::
Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV start_X,T_ARRAY_TV end_X)
{
    STATIC_ASSERT((IS_SAME<TV,typename T_ARRAY_TV::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY_TV>::value));
    if(triangles_per_group) for(int k=1;k<=leaves;k++){
        if(triangles_in_group(k).m){const VECTOR<int,3>& nodes=triangle_mesh.elements(triangles_in_group(k)(1));
            bounding_boxes(k)=RANGE<TV>::Bounding_Box(start_X.Subset(nodes));bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(end_X.Subset(nodes));}
        for(int i=2;i<=triangles_in_group(k).m;i++){const VECTOR<int,3>& nodes=triangle_mesh.elements(triangles_in_group(k)(i));
            bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(start_X.Subset(nodes));bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(end_X.Subset(nodes));}}
    else for(int k=1;k<=leaves;k++){const VECTOR<int,3>& nodes=triangle_mesh.elements(k);
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(start_X.Subset(nodes));bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(end_X.Subset(nodes));}
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class T> void TRIANGLE_HIERARCHY<T>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const FRAME<TV>& start_frame,const FRAME<TV>& end_frame)
{
    if(triangles_per_group) for(int k=1;k<=leaves;k++){
        if(triangles_in_group(k).m){const VECTOR<int,3>& nodes=triangle_mesh.elements(triangles_in_group(k)(1));
            bounding_boxes(k)=RANGE<TV>::Bounding_Box(start_frame*particles.X(nodes[1]),start_frame*particles.X(nodes[2]),start_frame*particles.X(nodes[3]));
            bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(end_frame*particles.X(nodes[1]),end_frame*particles.X(nodes[2]),end_frame*particles.X(nodes[3]));}
        for(int i=2;i<=triangles_in_group(k).m;i++){const VECTOR<int,3>& nodes=triangle_mesh.elements(triangles_in_group(k)(i));
            bounding_boxes(k)=RANGE<TV>::Bounding_Box(start_frame*particles.X(nodes[1]),start_frame*particles.X(nodes[2]),start_frame*particles.X(nodes[3]));
            bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(end_frame*particles.X(nodes[1]),end_frame*particles.X(nodes[2]),end_frame*particles.X(nodes[3]));}}
    else for(int k=1;k<=leaves;k++){const VECTOR<int,3>& nodes=triangle_mesh.elements(k);
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(start_frame*particles.X(nodes[1]),start_frame*particles.X(nodes[2]),start_frame*particles.X(nodes[3]));
        bounding_boxes(k).Enlarge_Nonempty_Box_To_Include_Points(end_frame*particles.X(nodes[1]),end_frame*particles.X(nodes[2]),end_frame*particles.X(nodes[3]));}
}
//#####################################################################
// Function Calculate_Bounding_Box_Radii
//#####################################################################
// at the boxes center, but may be a tighter bound on the triangle than the box
template<class T> void TRIANGLE_HIERARCHY<T>::
Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) 
{
    if(triangles_per_group) for(int k=1;k<=leaves;k++){
        VECTOR<T,3> center=bounding_boxes(k).Center();int node1,node2,node3;triangle_mesh.elements(k).Get(node1,node2,node3);
        radius(k)=0;
        for(int i=1;i<=triangles_in_group(i).m;i++)
            radius(k)=max(radius(k),(particles.X(node1)-center).Magnitude_Squared(),(particles.X(node2)-center).Magnitude_Squared(),(particles.X(node3)-center).Magnitude_Squared());
        radius(k)=sqrt(radius(k));}
    else for(int k=1;k<=leaves;k++){
        VECTOR<T,3> center=bounding_boxes(k).Center();int node1,node2,node3;triangle_mesh.elements(k).Get(node1,node2,node3);
        radius(k)=sqrt(max((particles.X(node1)-center).Magnitude_Squared(),(particles.X(node2)-center).Magnitude_Squared(),(particles.X(node3)-center).Magnitude_Squared()));}
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool TRIANGLE_HIERARCHY<T>::
Intersection(RAY<VECTOR<T,3> >& ray,const T thickness_over_two,const bool use_ray_bounding_box) const
{
    if(use_ray_bounding_box) ray.Compute_Bounding_Box();
    return Intersection(root,ray,thickness_over_two*2,thickness_over_two,use_ray_bounding_box);
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool TRIANGLE_HIERARCHY<T>::
Intersection(RAY<VECTOR<T,3> >& ray,VECTOR<T,3>& weights,const T thickness_over_two,const bool use_ray_bounding_box) const
{
    if(use_ray_bounding_box) ray.Compute_Bounding_Box();
    if(!Intersection(root,ray,thickness_over_two*2,thickness_over_two,use_ray_bounding_box)) return false;
    int i,j,k;triangle_mesh.elements(ray.aggregate_id).Get(i,j,k);
    TRIANGLE_3D<T>(particles.X(i),particles.X(j),particles.X(k)).Closest_Point(ray.Point(ray.t_max),weights);
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool TRIANGLE_HIERARCHY<T>::
Intersection(const int box,RAY<VECTOR<T,3> >& ray,const T thickness,const T thickness_over_two,const bool use_ray_bounding_box) const
{
    if(use_ray_bounding_box && !ray.semi_infinite) if(!box_hierarchy(box).Intersection(ray.bounding_box,thickness_over_two)) return false;
     
    if(box_hierarchy(box).Outside(ray.endpoint,thickness)){
        RAY<VECTOR<T,3> > ray_temp;ray.Save_Intersection_Information(ray_temp);
        if(!INTERSECTION::Lazy_Intersects(ray,box_hierarchy(box),thickness_over_two)) return false; // does not intersect the box or start inside it
        ray.Restore_Intersection_Information(ray_temp);}

    // if at the lowest level, check the triangle and return 1 or 0
    if(Leaf(box)){
        if(triangles_per_group){
            bool intersected=false;
            if(triangle_list){
                for(int i=1;i<=triangles_in_group(box).m;i++){
                    RAY<VECTOR<T,3> > ray_temp=ray;
                    if(INTERSECTION::Lazy_Intersects(ray_temp,(*triangle_list)(triangles_in_group(box)(i)).Bounding_Box(),thickness_over_two) && 
                       INTERSECTION::Intersects(ray,(*triangle_list)(triangles_in_group(box)(i)),thickness_over_two)){
                        ray.aggregate_id=triangles_in_group(box)(i);intersected=true;}}
                if(intersected && use_ray_bounding_box) ray.Compute_Bounding_Box();
                return intersected;}
            else{
                bool intersected=false;
                for(int i=1;i<=triangles_in_group(box).m;i++){
                    int node1,node2,node3;triangle_mesh.elements(triangles_in_group(box)(i)).Get(node1,node2,node3);
                    TRIANGLE_3D<T> triangle(particles.X(node1),particles.X(node2),particles.X(node3));
                    RAY<VECTOR<T,3> > ray_temp=ray;
                    if(INTERSECTION::Lazy_Intersects(ray_temp,triangle.Bounding_Box(),thickness_over_two) && INTERSECTION::Intersects(ray,triangle,thickness_over_two)){
                        ray.aggregate_id=triangles_in_group(box)(i);intersected=true;}}
                if(intersected && use_ray_bounding_box) ray.Compute_Bounding_Box();
                return intersected;}}
        else{
            if(triangle_list){
                if(INTERSECTION::Intersects(ray,(*triangle_list)(box),thickness_over_two)){
                    if(use_ray_bounding_box) ray.Compute_Bounding_Box();
                    ray.aggregate_id=box;return true;}
                else return false;}
            else{
                int node1,node2,node3;triangle_mesh.elements(box).Get(node1,node2,node3);
                TRIANGLE_3D<T> triangle(particles.X(node1),particles.X(node2),particles.X(node3));
                if(INTERSECTION::Intersects(ray,triangle,thickness_over_two)){
                    if(use_ray_bounding_box) ray.Compute_Bounding_Box();
                    ray.aggregate_id=box;return true;}
                else return false;}}}
    
    // if not at the lowest level, check the child boxes
    int box1,box2;children(box-leaves).Get(box1,box2);
    bool check1=Intersection(box1,ray,thickness,thickness_over_two),check2=Intersection(box2,ray,thickness,thickness_over_two);
    return check1 || check2;
}
//#####################################################################
template class TRIANGLE_HIERARCHY<float>;
template void TRIANGLE_HIERARCHY<float>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,ARRAY_VIEW<const TV>);
template void TRIANGLE_HIERARCHY<float>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,INDIRECT_ARRAY<ARRAY_VIEW<const TV> >);
template void TRIANGLE_HIERARCHY<float>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,ARRAY_VIEW<const TV>,ARRAY_VIEW<const TV>);
template void TRIANGLE_HIERARCHY<float>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,INDIRECT_ARRAY<ARRAY_VIEW<const TV> >,INDIRECT_ARRAY<ARRAY_VIEW<const TV> >);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_HIERARCHY<double>;
template void TRIANGLE_HIERARCHY<double>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,ARRAY_VIEW<const TV>);
template void TRIANGLE_HIERARCHY<double>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,INDIRECT_ARRAY<ARRAY_VIEW<const TV> >);
template void TRIANGLE_HIERARCHY<double>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,ARRAY_VIEW<const TV>,ARRAY_VIEW<const TV>);
template void TRIANGLE_HIERARCHY<double>::Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >&,INDIRECT_ARRAY<ARRAY_VIEW<const TV> >,INDIRECT_ARRAY<ARRAY_VIEW<const TV> >);
#endif
