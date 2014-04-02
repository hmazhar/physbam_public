//#####################################################################
// Copyright 2002-2008, Chris Allocco, Robert Bridson, Kevin Der, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Igor Neverov, Andrew Selle, Eftychios Sifakis, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGULATED_SURFACE
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SPHERE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_3D_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_CLEANSING.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_INSIDE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_INSIDE_USING_RAY_TEST.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_NORMAL.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_REFRESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_UPDATE_VERTEX_NORMALS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
using namespace PhysBAM;
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Triangulated_Surface(){
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<TRIANGULATED_SURFACE<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<TRIANGULATED_SURFACE<double> >();
#endif
    return true;
}
static bool registered=Register_Triangulated_Surface();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>::
TRIANGULATED_SURFACE(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,TRIANGLE_MESH>(triangle_mesh_input,particles_input),
    triangle_list(0),segment_lengths(0),hierarchy(0),particle_hierarchy(0),
    avoid_normal_interpolation_across_sharp_edges(false),normal_variance_threshold((T).1),
    vertex_normals(0),face_vertex_normals(0)
{
    PHYSBAM_ASSERT(registered);
    Use_Face_Normals();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>::
~TRIANGULATED_SURFACE()
{
    delete triangle_list;delete hierarchy;delete segment_lengths;delete particle_hierarchy;
    delete vertex_normals;delete face_vertex_normals;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,TRIANGLE_MESH>::Clean_Memory();
    delete triangle_list;triangle_list=0;delete hierarchy;hierarchy=0;
    delete segment_lengths;segment_lengths=0;delete particle_hierarchy;particle_hierarchy=0;delete vertex_normals;vertex_normals=0;
    delete face_vertex_normals;face_vertex_normals=0;
    Use_Face_Normals();
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Refresh_Auxiliary_Structures_Helper()
{
    if(triangle_list) Update_Triangle_List();
    if(hierarchy) Initialize_Hierarchy();
    if(vertex_normals || face_vertex_normals) TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Vertex_Normals(*this);
    if(segment_lengths) TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Segment_Lengths(*this);
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Hierarchy(const bool update_boxes,const int triangles_per_group) // creates and updates the boxes as well
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Hierarchy(*this,update_boxes,triangles_per_group);
}
//#####################################################################
// Function Initialize_Particle_Hierarchy
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Particle_Hierarchy(const INDIRECT_ARRAY<ARRAY_VIEW<TV> >& particle_subset_input,const bool update_boxes,const int particles_per_group) // creates and updates the boxes as well
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Particle_Hierarchy(*this,particle_subset_input,update_boxes,particles_per_group);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Rescale(const T scaling_x,const T scaling_y,const T scaling_z)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Rescale(*this,scaling_x,scaling_y,scaling_z);
}
//#####################################################################
// Function Update_Triangle_List
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Update_Triangle_List()
{
    Update_Triangle_List(particles.X);
}
//#####################################################################
// Function Update_Triangle_List
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Update_Triangle_List(ARRAY_VIEW<const TV> X)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Triangle_List(*this,X);
}
//#####################################################################
// Function Initialize_Torus_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Torus_Mesh_And_Particles(const int m,const int n,const T major_radius,const T minor_radius)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Torus_Mesh_And_Particles(*this,m,n,major_radius,minor_radius);
}
//#####################################################################
// Function Initialize_Cylinder_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Cylinder_Mesh_And_Particles(const int m,const int n,const T length,const T radius,const bool create_caps)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Cylinder_Mesh_And_Particles(*this,m,n,length,radius,create_caps);
}
//#####################################################################
// Function Initialize_Segment_Lengths
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Segment_Lengths()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Segment_Lengths(*this);
}
//#####################################################################
// Function Update_Vertex_Normals
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Update_Vertex_Normals()
{  
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Vertex_Normals(*this);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Normal(const TV& location,const int aggregate) const 
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Normal(*this,location,aggregate);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside(const TV& location,const T thickness_over_two) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Inside(*this,location,thickness_over_two);
}
//#####################################################################
// Function Inside_Relative_To_Triangle
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside_Relative_To_Triangle(const TV& location,const int triangle_index_for_ray_test,const T thickness_over_two) const
{
    PHYSBAM_ASSERT(triangle_list);
    RAY<TV> ray(SEGMENT_3D<T>(location,(*triangle_list)(triangle_index_for_ray_test).Center()));
    ray.t_max+=2*thickness_over_two;
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Inside_Using_Ray_Test(*this,ray,thickness_over_two); 
}
//#####################################################################
// Function Inside_Using_Ray_Test
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside_Using_Ray_Test(RAY<TV>& ray,const T thickness_over_two) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Inside_Using_Ray_Test(*this,ray,thickness_over_two); 
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Outside(const TV& location,const T thickness_over_two) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Outside(*this,location,thickness_over_two);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Inside_Any_Triangle
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside_Any_Triangle(const TV& location,int& triangle_id,const T thickness_over_two) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Inside_Any_Triangle(*this,location,triangle_id,thickness_over_two);
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Surface(const TV& location,const T max_depth,const T thickness_over_2,int* closest_triangle,T* distance) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Surface(*this,location,max_depth,thickness_over_2,closest_triangle,distance);
}
//#####################################################################
// Function Oriented_Surface
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Oriented_Surface(const TV& location,const TV& normal,const T max_depth,const T thickness_over_2,int* closest_triangle,T* distance) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Oriented_Surface(*this,location,normal,max_depth,thickness_over_2,closest_triangle,distance);
}
//#####################################################################
// Function Signed_Solid_Angle_Of_Triangle_Web
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Signed_Solid_Angle_Of_Triangle_Web(const TV& location,int web_root_node) const
{
    assert(mesh.incident_elements);
    T signed_solid_angle=0;
    for(int i=1;i<=(*mesh.incident_elements)(web_root_node).m;i++){
        int t=(*mesh.incident_elements)(web_root_node)(i);
        int node1,node2,node3;mesh.elements(t).Get(node1,node2,node3);
        TV x1=particles.X(node1),x2=particles.X(node2),x3=particles.X(node3);
        SPHERE<TV> sphere(location,1);
        // extend two of the triangle edges so their farther endpoints are on the unit sphere
        if(web_root_node == node1){RAY<TV> ray1(x1,x2-x1);INTERSECTION::Intersects(ray1,sphere,(T)0);x2=ray1.Point(ray1.t_max);RAY<TV> ray2(x1,x3-x1);INTERSECTION::Intersects(ray2,sphere,(T)0);x3=ray2.Point(ray2.t_max);} 
        else if(web_root_node == node2){RAY<TV> ray1(x2,x1-x2);INTERSECTION::Intersects(ray1,sphere,(T)0);x1=ray1.Point(ray1.t_max);RAY<TV> ray2(x2,x3-x2);INTERSECTION::Intersects(ray2,sphere,(T)0);x3=ray2.Point(ray2.t_max);} 
        else{RAY<TV> ray1(x3,x1-x3);INTERSECTION::Intersects(ray1,sphere,(T)0);x1=ray1.Point(ray1.t_max);RAY<TV> ray2(x3,x2-x3);INTERSECTION::Intersects(ray2,sphere,(T)0);x2=ray2.Point(ray2.t_max);}
        signed_solid_angle+=TRIANGLE_3D<T>(x1,x2,x3).Signed_Solid_Angle(location);} 
    return signed_solid_angle;
}
//#####################################################################
// Function Check_For_Self_Intersection
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Check_For_Self_Intersection(const T thickness_over_2,const bool update_bounding_boxes,ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs)
{  
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Check_For_Self_Intersection(*this,thickness_over_2,update_bounding_boxes,intersecting_segment_triangle_pairs);
}
//#####################################################################
// Function Find_First_Segment_Triangle_Intersection
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Find_First_Segment_Triangle_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const int max_coarsening_attempts,
    const bool update_bounding_boxes)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Find_First_Segment_Triangle_Intersection(*this,test_segment_mesh,X,thickness_over_2,max_coarsening_attempts,update_bounding_boxes);
}
//#####################################################################
// Function Segment_Triangle_Intersection
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Segment_Triangle_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const bool update_bounding_boxes,
    ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs)
{  
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Segment_Triangle_Intersection(*this,test_segment_mesh,X,thickness_over_2,update_bounding_boxes,intersecting_segment_triangle_pairs);
}
//#####################################################################
// Function Get_Triangles_Near_Edges
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Get_Triangles_Near_Edges(ARRAY<ARRAY<int> >& triangles_near_edges,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,
    const bool update_bounding_boxes)
{  
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Get_Triangles_Near_Edges(*this,triangles_near_edges,test_segment_mesh,X,thickness_over_2,update_bounding_boxes);
}
//#####################################################################
// Function Centroid_Of_Neighbors
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Centroid_Of_Neighbors(const int node) const
{
    assert(mesh.neighbor_nodes);
    TV target;
    int number_of_neighbors=(*mesh.neighbor_nodes)(node).m;
    for(int j=1;j<=number_of_neighbors;j++) target+=particles.X((*mesh.neighbor_nodes)(node)(j));
    if(number_of_neighbors != 0) target/=(T)number_of_neighbors;
    else target=particles.X(node); // if no neighbors, return the current node location
    return target;
}
//#####################################################################
// Function Calculate_Signed_Distance
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Calculate_Signed_Distance(const TV& location,T thickness) const
{
    TV closest_point=Surface(location); // this uses the slow (but accurate) method
    T distance=(closest_point-location).Magnitude();
    if(Inside(location,thickness)) distance*=-1;
    return distance;
}
//#####################################################################
// Function Linearly_Subdivide
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Linearly_Subdivide()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Linearly_Subdivide(*this);
}
//#####################################################################
// Function Loop_Subdivide
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Loop_Subdivide()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Loop_Subdivide(*this);
}
//#####################################################################
// Function Root_Three_Subdivide
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Root_Three_Subdivide()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Root_Three_Subdivide(*this);
}
//#####################################################################
// Funcion Total_Area
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Total_Area() const
{
    T area=0;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2),node3=mesh.elements(t)(3);
        area+=TRIANGLE_3D<T>::Area(particles.X(node1),particles.X(node2),particles.X(node3));}
    return area;
}
//#####################################################################
// Function Minimum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Angle(int* index) const
{
    TV u,v,w;T min_cosine=1;
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        T local_min=min(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u));
        if(local_min < min_cosine){min_cosine=local_min;if(index) *index=t;}}
    return acos(max(-min_cosine,(T)-1));
}
//#####################################################################
// Function Maximum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Angle(int* index) const
{
    TV u,v,w;T max_cosine=-1;
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        T local_max=max(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u));
        if(local_max > max_cosine){max_cosine=local_max;if(index) *index=t;}}
    return acos(min(-max_cosine,(T)1));
}
//#####################################################################
// Function Average_Minimum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Minimum_Angle() const
{
    TV u,v,w;T total_min_angle=0;
    for(int t=1; t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        total_min_angle+=acos(max(-min(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u)),(T)-1));}
    return total_min_angle/=mesh.elements.m;
}
//#####################################################################
// Function Average_Maximum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Maximum_Angle() const
{
    TV u,v,w;T total_max_angle=0;
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        total_max_angle+=acos(min(-max(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u)),(T)1));}
    return total_max_angle/=mesh.elements.m;
}
//#####################################################################
// Function Minimum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Edge_Length(int* index) const
{
    T min_edge_length_squared=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T local_min_squared=min((particles.X(i)-particles.X(j)).Magnitude_Squared(),(particles.X(j)-particles.X(k)).Magnitude_Squared(),(particles.X(k)-particles.X(i)).Magnitude_Squared());
        if(local_min_squared < min_edge_length_squared){min_edge_length_squared=local_min_squared;if(index) *index=t;}}
    return sqrt(min_edge_length_squared);
}
//#####################################################################
// Function Maximum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Edge_Length(int* index) const
{
    T max_edge_length_squared=0;
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T local_max_squared=max((particles.X(i)-particles.X(j)).Magnitude_Squared(),(particles.X(j)-particles.X(k)).Magnitude_Squared(),(particles.X(k)-particles.X(i)).Magnitude_Squared());
        if(local_max_squared > max_edge_length_squared){max_edge_length_squared=local_max_squared;if(index) *index=t;}}
    return sqrt(max_edge_length_squared);
}
//#####################################################################
// Function Average_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Edge_Length() const
{
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    T average_edge_length=0;
    for(int s=1;s<=mesh.segment_mesh->elements.m;s++){
        int i,j;mesh.segment_mesh->elements(s).Get(i,j);
        average_edge_length+=(particles.X(i)-particles.X(j)).Magnitude();}
    if(mesh.segment_mesh->elements.m) average_edge_length/=mesh.segment_mesh->elements.m;
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    return average_edge_length;
}
//#####################################################################
// Function Maximum_Aspect_Ratio
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Aspect_Ratio(int* index) const
{
    T max_aspect_ratio=0;
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T local_max=TRIANGLE_3D<T>::Aspect_Ratio(particles.X(i),particles.X(j),particles.X(k));
        if(local_max > max_aspect_ratio){max_aspect_ratio=local_max;if(index) *index=t;}}
    return max_aspect_ratio;
}
//#####################################################################
// Function Average_Aspect_Ratio
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Aspect_Ratio() const
{
    T total_aspect_ratio=0;
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        total_aspect_ratio+=TRIANGLE_3D<T>::Aspect_Ratio(particles.X(i),particles.X(j),particles.X(k));}
    return total_aspect_ratio/mesh.elements.m;
}
//#####################################################################
// Funcion Minimum_Area
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Area(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1,node2,node3;mesh.elements(t).Get(node1,node2,node3);
        T temp=TRIANGLE_3D<T>::Area(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Altitude
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Altitude(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2),node3=mesh.elements(t)(3);
        T temp=TRIANGLE_3D<T>::Minimum_Altitude(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Function Maximum_Magnitude_Phi
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Magnitude_Phi(const IMPLICIT_OBJECT<TV>& implicit_surface,int* index)
{
    T phi=0,max_phi=0;int k=0;
    for(int i=1;i<=particles.array_collection->Size();i++){phi=abs(implicit_surface(particles.X(i)));if(phi > max_phi){max_phi=phi;k=i;}}
    if(index)*index=k;return max_phi;
}
//#####################################################################
// Function Make_Orientations_Consistent_With_Implicit_Surface
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Make_Orientations_Consistent_With_Implicit_Surface(const IMPLICIT_OBJECT<TV>& implicit_surface)
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        TV centroid=(T)one_third*(particles.X(i)+particles.X(j)+particles.X(k));
        TV normal=TRIANGLE_3D<T>::Normal(particles.X(i),particles.X(j),particles.X(k));
        if(TV::Dot_Product(normal,implicit_surface.Normal(centroid))<0)mesh.elements(t).Set(i,k,j);}
}
//#####################################################################
// Function Close_Surface
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Close_Surface(const bool merge_coincident_vertices,const T merge_coincident_vertices_threshold,const bool fill_holes,const bool verbose)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Close_Surface(*this,merge_coincident_vertices,merge_coincident_vertices_threshold,fill_holes,verbose);
}
//#####################################################################
// Function Remove_Degenerate_Triangles
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Remove_Degenerate_Triangles(const T area_threshold)
{
    T area_threshold_squared=sqr(area_threshold);
    for(int t=mesh.elements.m;t>=1;t--){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        if(TRIANGLE_3D<T>::Area_Squared(particles.X(i),particles.X(j),particles.X(k))<area_threshold_squared) mesh.elements.Remove_Index_Lazy(t);}
    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Create_Compact_Copy
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* TRIANGULATED_SURFACE<T>::
Create_Compact_Copy() const
{
    ARRAY<int> new_to_old;
    HASHTABLE<int,int> old_to_new;
    TRIANGLE_MESH* triangle_mesh=new TRIANGLE_MESH;
    triangle_mesh->elements.Resize(mesh.elements.m);
    for(int i=1;i<=mesh.elements.m;i++){const VECTOR<int,3>& element=mesh.elements(i);
        for(int j=1;j<=3;j++){
            int& a=old_to_new.Get_Or_Insert(element(j));
            if(!a) a=new_to_old.Append(element(j));
            triangle_mesh->elements(i)(j)=a;}}

    GEOMETRY_PARTICLES<TV>* deformable_geometry_particle=new GEOMETRY_PARTICLES<TV>;
    deformable_geometry_particle->array_collection->Add_Elements(new_to_old.Size());
    deformable_geometry_particle->X.Prefix(deformable_geometry_particle->array_collection->Size())=particles.X.Subset(new_to_old);
    TRIANGULATED_SURFACE* triangulated_surface=new TRIANGULATED_SURFACE(*triangle_mesh,*deformable_geometry_particle);
    triangulated_surface->Update_Number_Nodes();
    return triangulated_surface;
}
//#####################################################################
template class TRIANGULATED_SURFACE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGULATED_SURFACE<double>;
#endif
