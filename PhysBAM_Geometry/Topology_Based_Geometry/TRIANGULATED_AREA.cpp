//#####################################################################
// Copyright 2003-2008, Christopher Allocco, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Neil Molino, Duc Nguyen, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_AREA_INSIDE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_AREA_PRUNE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_AREA_REFRESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_AREA_SPLIT.h>
using namespace PhysBAM;
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Triangulated_Area(){
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<TRIANGULATED_AREA<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<TRIANGULATED_AREA<double> >();
#endif
    return true;
}
static bool registered=Register_Triangulated_Area();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGULATED_AREA<T>::
TRIANGULATED_AREA(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<VECTOR<T,2> >& particles_input)
    :MESH_OBJECT<TV,TRIANGLE_MESH>(triangle_mesh_input,particles_input),segmented_curve(0),hierarchy(0),triangle_area_fractions(0),triangle_areas(0),
    nodal_areas(0)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TRIANGULATED_AREA<T>::
~TRIANGULATED_AREA()
{
    delete segmented_curve;delete hierarchy;delete triangle_area_fractions;delete triangle_areas;delete nodal_areas;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,TRIANGLE_MESH>::Clean_Memory();
    delete segmented_curve;segmented_curve=0;delete hierarchy;hierarchy=0;
    delete triangle_area_fractions;triangle_area_fractions=0;delete triangle_areas;triangle_areas=0;delete nodal_areas;nodal_areas=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Refresh_Auxiliary_Structures_Helper()
{
    if(hierarchy) Initialize_Hierarchy();
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Hierarchy(const bool update_boxes)
{
    delete hierarchy;hierarchy=new TRIANGLE_HIERARCHY_2D<T>(mesh,particles,update_boxes);
}
//#####################################################################
// Funcion Initialize_Square_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Square_Mesh_And_Particles(const GRID<TV>& grid,const bool reverse_triangles)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Square_Mesh_And_Particles(*this,grid,reverse_triangles);
}
//#####################################################################
// Funcion Initialize_Circle_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Circle_Mesh_And_Particles(const T outer_radius,const T inner_radius,const int num_radial,const int num_tangential)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Circle_Mesh_And_Particles(*this,outer_radius,inner_radius,num_radial,num_tangential);
}
//#####################################################################
// Funcion Initialize_Herringbone_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Herring_Bone_Mesh_And_Particles(const GRID<TV>& grid)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Herring_Bone_Mesh_And_Particles(*this,grid);
}
//#####################################################################
// Funcion Initialize_Equilateral_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Equilateral_Mesh_And_Particles(const GRID<TV>& grid)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Equilateral_Mesh_And_Particles(*this,grid);
}
//#####################################################################
// Function Inside
// note: return the first triangle that it is inside of (including boundary), otherwise returns 0
//#####################################################################
template<class T> int TRIANGULATED_AREA<T>::
Inside(const VECTOR<T,2>& location,const T thickness_over_two) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Inside(*this,location,thickness_over_two);
}
//#####################################################################
// Function Check_Signed_Area_And_Make_Consistent
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Check_Signed_Area_And_Make_Consistent(const int triangle,const bool verbose)
{
    VECTOR<int,3>& nodes=mesh.elements(triangle);
    if(TRIANGLE_2D<T>::Signed_Area(particles.X(nodes[1]),particles.X(nodes[2]),particles.X(nodes[3])) < 0){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"triangle number "<<triangle<<" is oriented improperly."<<std::endl;
#endif
        exchange(nodes[2],nodes[3]);}
}
//#####################################################################
// Function Check_Signed_Areas_And_Make_Consistent
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Check_Signed_Areas_And_Make_Consistent(const bool verbose)
{
    for(int t=1;t<=mesh.elements.m;t++) Check_Signed_Area_And_Make_Consistent(t,verbose);
}
//#####################################################################
// Function Initialize_Segmented_Curve
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Segmented_Curve()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Segmented_Curve(*this);
}
//#####################################################################
// Function Centroid
//#####################################################################
template<class T> VECTOR<T,2> TRIANGULATED_AREA<T>::
Centroid(const int triangle) const
{
    int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return (T)one_third*(particles.X(i)+particles.X(j)+particles.X(k));
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Rescale(const T scaling_factor)
{
    Rescale(scaling_factor,scaling_factor);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Rescale(const T scaling_x,const T scaling_y)
{
    for(int k=1;k<=particles.array_collection->Size();k++) particles.X(k)*=VECTOR<T,2>(scaling_x,scaling_y);
    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Area(const int triangle) const
{
    int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return TRIANGLE_2D<T>::Area(particles.X(i),particles.X(j),particles.X(k));
}
//#####################################################################
// Function Signed_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Signed_Area(const int triangle) const
{
    int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return TRIANGLE_2D<T>::Signed_Area(particles.X(i),particles.X(j),particles.X(k));
}
//#####################################################################
// Funcion Minimum_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Area(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2),node3=mesh.elements(t)(3);
        T temp=TRIANGLE_2D<T>::Area(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Signed_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Signed_Area(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2),node3=mesh.elements(t)(3);
        T temp=TRIANGLE_2D<T>::Signed_Area(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Total_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Total_Area() const
{
    T area=0;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2),node3=mesh.elements(t)(3);
        area+=TRIANGLE_2D<T>::Area(particles.X(node1),particles.X(node2),particles.X(node3));}
    return area;
}
//#####################################################################
// Funcion Minimum_Altitude
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Altitude(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2),node3=mesh.elements(t)(3);
        T temp=TRIANGLE_2D<T>::Minimum_Altitude(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Edge_Length(int* index) const
{
    int t_save=0;T minimum_squared=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){int i,j,k;mesh.elements(t).Get(i,j,k);
        T temp=TRIANGLE_2D<T>::Minimum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k));if(temp < minimum_squared){minimum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(minimum_squared);
}
//#####################################################################
// Funcion Maximum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Maximum_Edge_Length(int* index) const
{
    int t_save=0;T maximum_squared=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){int i,j,k;mesh.elements(t).Get(i,j,k);
        T temp=TRIANGLE_2D<T>::Maximum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k));if(temp < maximum_squared){maximum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(maximum_squared);
}
//#####################################################################
// Function Inverted_Triangles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Inverted_Triangles(ARRAY<int>& inverted_triangles) const
{
    inverted_triangles.Resize(0);
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2),node3=mesh.elements(t)(3);
        if(TRIANGLE_2D<T>::Signed_Area(particles.X(node1),particles.X(node2),particles.X(node3)) < 0) inverted_triangles.Append(t);}
}
//#####################################################################
// Function Volume_Incident_On_A_Particle
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Area_Incident_On_A_Particle(const int particle_index)
{
    int incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    T total_incident_area=0;
    for(int t=1;t<=(*mesh.incident_elements)(particle_index).m;t++){
        int i,j,k;mesh.elements((*mesh.incident_elements)(particle_index)(t)).Get(i,j,k);
        total_incident_area+=TRIANGLE_2D<T>::Area(particles.X(i),particles.X(j),particles.X(k));}
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return total_incident_area;
}
//#####################################################################
// Function Split_Node
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements
template<class T> int TRIANGULATED_AREA<T>::
Split_Node(const int particle_index,const VECTOR<T,2>& normal)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Split_Node(*this,particle_index,normal);
}
//#####################################################################
// Function Split_Connected_Component
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements & adjacent_elements
template<class T> int TRIANGULATED_AREA<T>::
Split_Connected_Component(const int node)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Split_Connected_Component(*this,node);
}
//#####################################################################
// Function Discard_Triangles_Outside_Implicit_Curve
//#####################################################################
// uses Whitney-like criterion to discard only those triangles that are for sure outside levelset (assuming accurate signed distance)
template<class T> void TRIANGULATED_AREA<T>::
Discard_Triangles_Outside_Implicit_Curve(IMPLICIT_OBJECT<TV>& implicit_curve)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Discard_Triangles_Outside_Implicit_Curve(*this,implicit_curve);
}
//#####################################################################
// Function Initialize_Triangle_Area_Fractions_From_Voronoi_Regions
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Triangle_Area_Fractions_From_Voronoi_Regions()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Triangle_Area_Fractions_From_Voronoi_Regions(*this);
}
//#####################################################################
// Function Compute_Triangle_Areas
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Compute_Triangle_Areas()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Compute_Triangle_Areas(*this);
}
//#####################################################################
// Function Compute_Nodal_Areas
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Compute_Nodal_Areas(bool save_triangle_areas)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Compute_Nodal_Areas(*this,save_triangle_areas);
}
//#####################################################################
// Function Triangle_In_Direction_Uninverted
//#####################################################################
template<class T> int TRIANGULATED_AREA<T>::
Triangle_In_Direction_Uninverted(const int node,const VECTOR<T,2>& direction) const
{
    assert(mesh.topologically_sorted_neighbor_nodes && mesh.topologically_sorted_incident_elements);
    ARRAY<int> &neighbor_nodes=(*mesh.topologically_sorted_neighbor_nodes)(node),&incident_elements=(*mesh.topologically_sorted_incident_elements)(node);
    VECTOR<T,2> xnode=particles.X(node);int t=1;
    while(t<neighbor_nodes.m){ // find right edge in right half-space of direction
        if(VECTOR<T,2>::Cross_Product(particles.X(neighbor_nodes(t))-xnode,direction).x<0){t++;continue;}
        while(t<neighbor_nodes.m){ // find left edge in left half-space of direction
            if(VECTOR<T,2>::Cross_Product(direction,particles.X(neighbor_nodes(t+1))-xnode).x<0){t++;continue;}
            return t;}
        return t<=incident_elements.m?incident_elements(t):0;}
    return t<=incident_elements.m?incident_elements(t):0;
}
//#####################################################################
// Function Triangle_Walk_Uninverted
//#####################################################################
template<class T> int TRIANGULATED_AREA<T>::
Triangle_Walk_Uninverted(const int start_node,const VECTOR<T,2>& dX) const
{
    int t=Triangle_In_Direction_Uninverted(start_node,dX);if(!t)return 0;
    VECTOR<T,2> start=particles.X(start_node),goal=start+dX;
    int e1,e2;mesh.Other_Two_Nodes(start_node,t,e1,e2);
    if(VECTOR<T,2>::Cross_Product(particles.X(e2)-particles.X(e1),goal-particles.X(e1)).x>=0)return t;
    for(;;){
        t=mesh.Adjacent_Triangle(t,e1,e2);if(!t)return 0;
        int e3=mesh.Other_Node(e1,e2,t);
        VECTOR<T,2> w=TRIANGLE_2D<T>::First_Two_Barycentric_Coordinates(goal,particles.X(e2),particles.X(e1),particles.X(e3));
        if(w.x>=0){if(w.y>=0)return t;else e1=e3;}
        else if(w.y>=0 || VECTOR<T,2>::Cross_Product(dX,particles.X(e3)-start).x>0) e2=e3;
        else e1=e3;}
}
//#####################################################################
// Function Triangle_Walk_Uninverted
//#####################################################################
template<class T> bool TRIANGULATED_AREA<T>::
Fix_Pair_For_Delaunay(const int triangle1,const int triangle2)
{
    VECTOR<int,3> &t1=mesh.elements(triangle1),&t2=mesh.elements(triangle2);
    int a1i=0;for(int j=1;j<=3;j++) if(!t2.Contains(t1(j))) a1i=j;
    int a2i=0;for(int j=1;j<=3;j++) if(!t1.Contains(t2(j))) a2i=j;
    int a1=t1(a1i),a2=t2(a2i),b1=t1(a1i>2?a1i-2:a1i+1),b2=t1(a1i>1?a1i-1:a1i+2);

    assert(TRIANGLE_2D<T>(particles.X.Subset(t1)).Check_Orientation() && TRIANGLE_2D<T>(particles.X.Subset(t2)).Check_Orientation());
    assert(b2==t2(a2i>2?a2i-2:a2i+1) && b1==t2(a2i>1?a2i-1:a2i+2));

    if(TRIANGLE_2D<T>::Check_Delaunay_Criterion(particles.X(a1),particles.X(b1),particles.X(b2),particles.X(a2))) return false;

    t1=VECTOR<int,3>(a1,b1,a2);t2=VECTOR<int,3>(a2,b2,a1);
    assert(TRIANGLE_2D<T>(particles.X.Subset(t1)).Check_Orientation() && TRIANGLE_2D<T>(particles.X.Subset(t2)).Check_Orientation());
    return true;
}
//#####################################################################
template class TRIANGULATED_AREA<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGULATED_AREA<double>;
#endif
