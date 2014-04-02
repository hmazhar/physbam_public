//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ron Fedkiw, Eilene Hao, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Mike Rodgers, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRALIZED_VOLUME
//#####################################################################
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TETRAHEDRALIZED_VOLUME_CENTROID.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TETRAHEDRALIZED_VOLUME_PRUNE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TETRAHEDRALIZED_VOLUME_REFRESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TETRAHEDRALIZED_VOLUME_SPLIT.h>
using namespace PhysBAM;
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Tetrahedralized_Volume(){
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<TETRAHEDRALIZED_VOLUME<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<TETRAHEDRALIZED_VOLUME<double> >();
#endif
    return true;
}
static bool registered=Register_Tetrahedralized_Volume();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRALIZED_VOLUME<T>::
TETRAHEDRALIZED_VOLUME(TETRAHEDRON_MESH& tetrahedron_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,TETRAHEDRON_MESH>(tetrahedron_mesh_input,particles_input),tetrahedron_list(0),
    triangulated_surface(0),hierarchy(0),tetrahedron_volumes(0),nodal_volumes(0)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TETRAHEDRALIZED_VOLUME<T>::
~TETRAHEDRALIZED_VOLUME()
{
    delete tetrahedron_list;delete triangulated_surface;delete hierarchy;delete tetrahedron_volumes;delete nodal_volumes;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,TETRAHEDRON_MESH>::Clean_Memory();
    delete tetrahedron_list;tetrahedron_list=0;delete triangulated_surface;triangulated_surface=0;
    delete hierarchy;hierarchy=0;delete tetrahedron_volumes;tetrahedron_volumes=0;delete nodal_volumes;nodal_volumes=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Refresh_Auxiliary_Structures_Helper()
{
    if(tetrahedron_list) Update_Tetrahedron_List();if(hierarchy) Initialize_Hierarchy();
    if(triangulated_surface) Initialize_Triangulated_Surface();if(tetrahedron_volumes) Compute_Tetrahedron_Volumes();
    if(nodal_volumes) Compute_Nodal_Volumes();
}
//#####################################################################
// Function Update_Tetrahedron_List
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Update_Tetrahedron_List()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Tetrahedron_List(*this);
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Hierarchy(const bool update_boxes) // creates and updates the boxes as well
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Hierarchy(*this,update_boxes);
}
//#####################################################################
// Funcion Initialize_Octahedron_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Octahedron_Mesh_And_Particles(const GRID<TV>& grid)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Octahedron_Mesh_And_Particles(*this,grid);
}
//#####################################################################
// Funcion Initialize_Cube_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Cube_Mesh_And_Particles(const GRID<TV>& grid)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Cube_Mesh_And_Particles(*this,grid);
}
//#####################################################################
// Function Initialize_Prismatic_Cube_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Prismatic_Cube_Mesh_And_Particles(const GRID<TV>& grid)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Prismatic_Cube_Mesh_And_Particles(*this,grid);
}
//#####################################################################
// Function Check_Signed_Volumes_And_Make_Consistent
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Check_Signed_Volumes_And_Make_Consistent(bool verbose)
{
    for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TV x1(particles.X(i)),x2(particles.X(j)),x3(particles.X(k)),x4(particles.X(l));
        T sign_of_volume=TV::Dot_Product(TV::Cross_Product(x2-x1,x3-x1),x4-x1); // left out division by 6
        if(sign_of_volume < 0){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose) LOG::cout<<"tetrahedron number "<<t<<" is oriented improperly."<<std::endl;
#endif
            exchange(mesh.elements(t)(3),mesh.elements(t)(4));}}
    if(tetrahedron_list) Update_Tetrahedron_List();
}
//#####################################################################
// Function Initialize_Triangulated_Surface
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Triangulated_Surface()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Triangulated_Surface(*this);
}
//#####################################################################
// Funcion Minimum_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Volume(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Volume();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Signed_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Signed_Volume(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Signed_Volume();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Total_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Total_Volume() const
{
    T volume=0;
    for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        volume+=TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));}
    return volume;
}
//#####################################################################
// Funcion Minimum_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Angle(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Minimum_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Minimum_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Maximum_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Angle(int* index) const
{
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Maximum_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Maximum_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Minimum_Altitude
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Altitude(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Minimum_Altitude();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Minimum_Altitude();if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Edge_Length(int* index) const
{
    int t_save=0;T minimum_squared=FLT_MAX;
    for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Minimum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp < minimum_squared){minimum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(minimum_squared);
}
//#####################################################################
// Funcion Maximum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Edge_Length(int* index) const
{
    int t_save=0;T maximum_squared=0;
    for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Maximum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp > maximum_squared){maximum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(maximum_squared);
}
//#####################################################################
// Funcion Maximum_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Aspect_Ratio(int* index) const
{
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Maximum_Interior_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Interior_Aspect_Ratio(int* index)
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m == 4){
        T temp=(*tetrahedron_list)(t).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}}
    else for(int t=1;t<=mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m == 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Maximum_Boundary_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Boundary_Aspect_Ratio(int* index)
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m != 4){
        T temp=(*tetrahedron_list)(t).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}}
    else for(int t=1;t<=mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m != 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Average_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Average_Aspect_Ratio()
{
    int total=0;T sum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){total++;sum+=(*tetrahedron_list)(t).Aspect_Ratio();}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        total++;sum+=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();}
    if(total != 0) sum/=total;
    return sum;
}
//#####################################################################
// Funcion Average_Interior_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Average_Interior_Aspect_Ratio()
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int total=0;T sum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m == 4){total++;sum+=(*tetrahedron_list)(t).Aspect_Ratio();}}
    else for(int t=1;t<=mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m == 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        total++;sum+=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(total != 0) sum/=total;
    return sum;
}
//#####################################################################
// Funcion Average_Boundary_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Average_Boundary_Aspect_Ratio()
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int total=0;T sum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m != 4){total++;sum+=(*tetrahedron_list)(t).Aspect_Ratio();}}
    else for(int t=1;t<=mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m != 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        total++;sum+=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(total != 0) sum/=total;
    return sum;
}
//#####################################################################
// Funcion Minimum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Dihedral_Angle(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Minimum_Dihedral_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Minimum_Dihedral_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Maximum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Dihedral_Angle(int* index) const
{
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=1;t<=mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Maximum_Dihedral_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    else for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Maximum_Dihedral_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Maximum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Edge_Length(int* index)
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    int s_save=0;T maximum=0;
    for(int s=1;s<=mesh.segment_mesh->elements.m;s++){int i,j;mesh.segment_mesh->elements(s).Get(i,j);
        T temp=(particles.X(i)-particles.X(j)).Magnitude();if(temp > maximum){maximum=temp;s_save=s;}}
    if(index) *index=s_save;
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    return maximum;
}
//#####################################################################
// Funcion Minimum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Edge_Length(int* index)
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    int s_save=0;T minimum=FLT_MAX;
    for(int s=1;s<=mesh.segment_mesh->elements.m;s++){int i,j;mesh.segment_mesh->elements(s).Get(i,j);
        T temp=(particles.X(i)-particles.X(j)).Magnitude();if(temp < minimum){minimum=temp;s_save=s;}}
    if(index) *index=s_save;
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    return minimum;
}
//#####################################################################
// Function Interior_Laplacian_Smoothing
//#####################################################################
// one step of mesh adjustment using Laplacian smoothing
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Advance_Interior_Laplacian_Smoothing()
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool neighbor_nodes_defined=mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined) mesh.Initialize_Neighbor_Nodes();
    bool node_on_boundary_defined=mesh.node_on_boundary!=0;if(!node_on_boundary_defined) mesh.Initialize_Node_On_Boundary();
    // compute the centroid of the neighbors - if on the boundary, just use boundary neighbors
    for(int i=1;i<=mesh.neighbor_nodes->m;i++){
        int number=0;TV target;
        for(int j=1;j<=(*mesh.neighbor_nodes)(i).m;j++){
            int node=(*mesh.neighbor_nodes)(i)(j);
            if(!(*mesh.node_on_boundary)(i) || (*mesh.node_on_boundary)(node)){number++;target+=particles.X(node);}}
        if(number != 0){target/=(T)number;particles.X(i)=target;}}
    if(!neighbor_nodes_defined){delete mesh.neighbor_nodes;mesh.neighbor_nodes=0;}
    if(!node_on_boundary_defined){delete mesh.node_on_boundary;mesh.node_on_boundary=0;}
}
//#####################################################################
// Function Centroid
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRALIZED_VOLUME<T>::
Centroid(const int tetrahedron) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return (T).25*(particles.X(i)+particles.X(j)+particles.X(k)+particles.X(l));
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Rescale(const T scaling_factor)
{
    Rescale(scaling_factor,scaling_factor,scaling_factor);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Rescale(const T scaling_x,const T scaling_y,const T scaling_z)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Rescale(*this,scaling_x,scaling_y,scaling_z);
}
//#####################################################################
// Function Signed_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Signed_Volume(const int tetrahedron) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Volume(const int tetrahedron) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
}
//#####################################################################
// Function Centroid_Of_Neighbors
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRALIZED_VOLUME<T>::
Centroid_Of_Neighbors(const int node) const
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Centroid_Of_Neighbors(*this,node);
}
//#####################################################################
// Function Completely_Inside_Box
//#####################################################################
template<class T> bool TETRAHEDRALIZED_VOLUME<T>::
Completely_Inside_Box(const int tetrahedron,const RANGE<TV>& box) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return (box.Lazy_Inside(particles.X(i)) && box.Lazy_Inside(particles.X(j)) && box.Lazy_Inside(particles.X(k)) && box.Lazy_Inside(particles.X(l)));
}
//#####################################################################
// Function Discard_Spikes_From_Adjacent_Elements
//#####################################################################
// throws out all tetrahedrons with only one neighbor (i.e. a spike on the boundary)
// returns index of first discarded
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Spikes_From_Adjacent_Elements(ARRAY<int>* deletion_list)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Discard_Spikes_From_Adjacent_Elements(*this,deletion_list);
}
//#####################################################################
// Function Interior_Edges_With_Boundary_Nodes
//#####################################################################
// throws out all tetrahedrons with all four nodes on boundary
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Interior_Edges_With_Boundary_Nodes(ARRAY<VECTOR<int,2> >* deletion_list)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Interior_Edges_With_Boundary_Nodes(*this,deletion_list);
}
//#####################################################################
// Function Discard_Spikes
//#####################################################################
// throws out all tetrahedrons with all four nodes on boundary
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Spikes(ARRAY<int>* deletion_list)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Discard_Spikes(*this,deletion_list);
}
//#####################################################################
// Function Inverted_Tetrahedrons
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Inverted_Tetrahedrons(ARRAY<int>& inverted_tetrahedrons) const
{
    inverted_tetrahedrons.Resize(0);
    for(int t=1;t<=mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TV x1(particles.X(i)),x2(particles.X(j)),x3(particles.X(k)),x4(particles.X(l));
        T sign_of_volume=TV::Dot_Product(TV::Cross_Product(x2-x1,x3-x1),x4-x1); // left out division by 6
        if(sign_of_volume < 0) inverted_tetrahedrons.Append(t);}
}
//#####################################################################
// Function Inside
//#####################################################################
// note: return the first tetrahedron that it is inside of (including boundary), otherwise returns 0
template<class T> int TETRAHEDRALIZED_VOLUME<T>::
Inside(const TV& location,const T thickness_over_two) const
{
    if(!hierarchy) PHYSBAM_FATAL_ERROR();
    if(hierarchy->box_hierarchy(hierarchy->root).Outside(location,thickness_over_two)) return 0;
        ARRAY<int> tetrahedrons_to_check;hierarchy->Intersection_List(location,tetrahedrons_to_check,thickness_over_two);
    if(tetrahedron_list) for(int p=1;p<=tetrahedrons_to_check.m;p++){
        int t=tetrahedrons_to_check(p);if(!(*tetrahedron_list)(t).Outside(location,thickness_over_two)) return t;}
    else for(int p=1;p<=tetrahedrons_to_check.m;p++){
        int t=tetrahedrons_to_check(p);int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TETRAHEDRON<T> tetrahedron_to_check(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
        if(!tetrahedron_to_check.Outside(location,thickness_over_two)) return t;}
    return 0;
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface
//#####################################################################
// uses Whitney-like criterion to discard only those tets that are for sure outside levelset (assuming accurate signed distance)
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Tetrahedrons_Outside_Implicit_Surface(IMPLICIT_OBJECT<TV>& implicit_surface)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Discard_Tetrahedrons_Outside_Implicit_Surface(*this,implicit_surface);
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive
//#####################################################################
// uses Whitney-like criterion to discard only those tets that are not for sure inside levelset (assuming accurate signed distance)
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(IMPLICIT_OBJECT<TV>& implicit_surface)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(*this,implicit_surface);
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive
//#####################################################################
// same as above Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive but restricted to tets fully contained inside one of the supplied boxes
// TODO: Need to decide what constitutes good criteria for tet to be contained in a box: for now demanding full containment, i.e. all four particles
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(IMPLICIT_OBJECT<TV>& implicit_surface,const ARRAY<RANGE<TV> >& bounding_boxes)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(*this,implicit_surface,bounding_boxes);
}
//#####################################################################
// Function Maximum_Magnitude_Phi_On_Boundary
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Magnitude_Phi_On_Boundary(const IMPLICIT_OBJECT<TV>& implicit_surface,int* index)
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool node_on_boundary_defined=mesh.node_on_boundary!=0;if(!node_on_boundary_defined) mesh.Initialize_Node_On_Boundary();
    T phi=0,max_phi=0;int p_save=0;
    for(int p=1;p<=particles.array_collection->Size();p++) if((*mesh.node_on_boundary)(p)){phi=abs(implicit_surface(particles.X(p)));if(phi > max_phi){max_phi=phi;p_save=p;}}
    if(index) *index=p_save;
    if(!node_on_boundary_defined){delete mesh.node_on_boundary;mesh.node_on_boundary=0;}
    return max_phi;
}
//#####################################################################
// Function Volume_Incident_On_A_Particle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Volume_Incident_On_A_Particle(const int particle_index)
{
    if(mesh.number_nodes!=particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    bool incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    T total_incident_volume=0;
    for(int t=1;t<=(*mesh.incident_elements)(particle_index).m;t++){int i,j,k,l;mesh.elements((*mesh.incident_elements)(particle_index)(t)).Get(i,j,k,l);
        total_incident_volume+=TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));}
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return total_incident_volume;
}
//#####################################################################
// Function Split_Along_Fracture_Plane
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Split_Along_Fracture_Plane(const PLANE<T>& plane,ARRAY<int>& particle_replicated)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Split_Along_Fracture_Plane(*this,plane,particle_replicated);
}
//#####################################################################
// Function Split_Node
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements
template<class T> int TETRAHEDRALIZED_VOLUME<T>::
Split_Node(const int particle_index,const TV& normal)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Split_Node(*this,particle_index,normal);
}
//#####################################################################
// Function Split_Connected_Component
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements & adjacent_elements
template<class T> int TETRAHEDRALIZED_VOLUME<T>::
Split_Connected_Component(const int node)
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Split_Connected_Component(*this,node);
}
//#####################################################################
// Function Compute_Tetrahedron_Volumes
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Compute_Tetrahedron_Volumes()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Compute_Tetrahedron_Volumes(*this);
}
//#####################################################################
// Function Compute_Nodal_Volumes
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Compute_Nodal_Volumes(bool save_tetrahedron_volumes)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Compute_Nodal_Volumes(*this,save_tetrahedron_volumes);
}
//#####################################################################
template class TETRAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TETRAHEDRALIZED_VOLUME<double>;
#endif
