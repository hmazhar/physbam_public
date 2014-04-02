//#####################################################################
// Copyright 2003-2009, Ron Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Duc Nguyen, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENTED_CURVE
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/SEGMENTED_CURVE_REFRESH.h>
namespace PhysBAM{
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Segmented_Curve(){
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<SEGMENTED_CURVE<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<SEGMENTED_CURVE<VECTOR<float,1> > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<SEGMENTED_CURVE<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<SEGMENTED_CURVE<VECTOR<double,1> > >();
#endif
    return true;
}
static bool registered=Register_Segmented_Curve();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SEGMENTED_CURVE<TV>::
SEGMENTED_CURVE(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,SEGMENT_MESH>(segment_mesh_input,particles_input),hierarchy(0),segment_list(0),point_simplices_1d(0)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SEGMENTED_CURVE<TV>::
~SEGMENTED_CURVE()
{
    delete hierarchy;delete segment_list;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Clean_Memory()
{
    MESH_OBJECT<TV,SEGMENT_MESH>::Clean_Memory();
    delete hierarchy;hierarchy=0;delete segment_list;segment_list=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Refresh_Auxiliary_Structures_Helper()
{
    if(segment_list) Update_Segment_List();
    if(hierarchy) Initialize_Hierarchy();
}
//#####################################################################
// Function Update_Segment_List
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Update_Segment_List() // updates the segments assuming the particle positions are already updated
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Segment_List(*this);
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::  
Initialize_Hierarchy(const bool update_boxes) // creates and updates the boxes as well
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Hierarchy(*this,update_boxes);
}
//#####################################################################
// Function Initialize_Straight_Mesh_And_Particles
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::  
Initialize_Straight_Mesh_And_Particles(const GRID<VECTOR<T,1> >& grid)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Straight_Mesh_And_Particles(*this,grid);
}
//#####################################################################
// Function Initialize_Circle_Mesh_And_Particles
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::  
Initialize_Circle_Mesh_And_Particles(const int m,const T radius)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Circle_Mesh_And_Particles(*this,m,radius);
}
template<> void SEGMENTED_CURVE<VECTOR<double,1> >::  
Initialize_Circle_Mesh_And_Particles(const int m,const T radius)
{PHYSBAM_NOT_IMPLEMENTED();}
template<> void SEGMENTED_CURVE<VECTOR<float,1> >::  
Initialize_Circle_Mesh_And_Particles(const int m,const T radius)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Closest_Point_On_Curve
//#####################################################################
template<class TV> TV SEGMENTED_CURVE<TV>::  
Closest_Point_On_Curve(const TV& location,T thickness_over_two,int* closest_segment,T* distance) const
{
    T min_distance_squared=FLT_MAX;TV point;
    for(int i=1;i<=mesh.elements.m;i++){
        T_SEGMENT segment=segment_list?(*segment_list)(i):Get_Element(i);
        TV new_point=segment.Closest_Point_On_Segment(location);
        T distance_squared=(new_point-location).Magnitude_Squared();
        if(distance_squared < min_distance_squared){
            min_distance_squared=distance_squared;point=new_point;
            if(closest_segment) *closest_segment=i;}}
    if(distance) *distance=sqrt(min_distance_squared);
    return point;
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Rescale(const TV& scaling)
{
    for(int k=1;k<=particles.array_collection->Size();k++) particles.X(k)*=scaling;
    if(segment_list) Update_Segment_List();if(bounding_box) Update_Bounding_Box();
}
//#####################################################################
// Function Average_Edge_Length
//#####################################################################
template<class TV> typename TV::SCALAR SEGMENTED_CURVE<TV>::
Average_Edge_Length() const
{
    T average_edge_length=0;
    for(int s=1;s<=mesh.elements.m;s++){
        int i,j;mesh.elements(s).Get(i,j);
        average_edge_length+=(particles.X(i)-particles.X(j)).Magnitude();}
    if(mesh.elements.m) average_edge_length/=mesh.elements.m;
    return average_edge_length;
}
//#####################################################################
// Function Total_Length
//#####################################################################
template<class TV> typename TV::SCALAR SEGMENTED_CURVE<TV>::
Total_Length() const
{
    T length=0;
    for(int t=1;t<=mesh.elements.m;t++){
        int node1=mesh.elements(t)(1),node2=mesh.elements(t)(2);
        length+=(particles.X(node1)-particles.X(node2)).Magnitude();}
    return length;
}
template<class T> void Initialize_Boundary_Object_Helper(SEGMENTED_CURVE<VECTOR<T,1> >& sc)
{
    delete sc.point_simplices_1d;if(!sc.mesh.boundary_mesh) sc.mesh.Initialize_Boundary_Mesh();
    sc.point_simplices_1d=new POINT_SIMPLICES_1D<T>(*sc.mesh.boundary_mesh,sc.particles);
}
template<class T> void Initialize_Boundary_Object_Helper(SEGMENTED_CURVE<VECTOR<T,2> >& sc){PHYSBAM_NOT_IMPLEMENTED();}
template<class T> void Initialize_Boundary_Object_Helper(SEGMENTED_CURVE<VECTOR<T,3> >& sc){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Get_Boundary_Object
//#####################################################################
template<class TV> POINT_SIMPLICES_1D<typename TV::SCALAR>& SEGMENTED_CURVE<TV>::
Get_Boundary_Object()
{
    Initialize_Boundary_Object_Helper(*this);
    return *point_simplices_1d;
}
//#####################################################################
template class SEGMENTED_CURVE<VECTOR<float,1> >;
template class SEGMENTED_CURVE<VECTOR<float,2> >;
template class SEGMENTED_CURVE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEGMENTED_CURVE<VECTOR<double,1> >;
template class SEGMENTED_CURVE<VECTOR<double,2> >;
template class SEGMENTED_CURVE<VECTOR<double,3> >;
#endif
}
