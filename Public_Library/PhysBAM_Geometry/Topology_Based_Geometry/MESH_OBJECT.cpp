//#####################################################################
// Copyright 2004-2009, Geoffrey Irving, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/FACTORIAL.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/MESH_OBJECT_APPEND.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/MESH_OBJECT_PRUNE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/MESH_OBJECT_REFRESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/MESH_OBJECT_UNION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T_MESH> MESH_OBJECT<TV,T_MESH>::
MESH_OBJECT(T_MESH& mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :mesh(mesh_input),particles(particles_input),bounding_box(0),particle_partition(0),need_destroy_mesh(false),need_destroy_particles(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T_MESH> MESH_OBJECT<TV,T_MESH>::
~MESH_OBJECT()
{
    if(need_destroy_mesh) delete &mesh;if(need_destroy_particles) delete &particles;
    delete bounding_box;
    delete particle_partition;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Clean_Memory()
{
    delete bounding_box;bounding_box=0;
    delete particle_partition;particle_partition=0;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,class T_MESH> typename MESH_TO_OBJECT<TV,T_MESH>::TYPE* MESH_OBJECT<TV,T_MESH>::
Create()
{
    T_DERIVED_OBJECT* object=new T_DERIVED_OBJECT(*(new T_MESH),*(new GEOMETRY_PARTICLES<TV>));
    object->need_destroy_mesh=object->need_destroy_particles=true;return object;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,class T_MESH> typename MESH_TO_OBJECT<TV,T_MESH>::TYPE* MESH_OBJECT<TV,T_MESH>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    T_DERIVED_OBJECT* object=new T_DERIVED_OBJECT(*(new T_MESH),particles);
    object->need_destroy_mesh=true;return object;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Refresh_Auxiliary_Structures()
{
    mesh.Refresh_Auxiliary_Structures();
    Refresh_Auxiliary_Structures_Helper();
    if(bounding_box) Update_Bounding_Box(); // this can depend on the hierarchy which by derived helper
    if(particle_partition) Initialize_Particle_Partition(particle_partition->grid.counts);
}
//#####################################################################
// Function Update_Number_Nodes
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Update_Number_Nodes()
{
    mesh.Set_Number_Nodes(particles.array_collection->Size());
}
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Update_Bounding_Box()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Bounding_Box(*this);
}
//#####################################################################
// Function Initialize_Particle_Partition
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Initialize_Particle_Partition(const VECTOR<int,TV::m>& counts)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Particle_Partition(*this,counts); 
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,class T_MESH> STRUCTURE<TV>* MESH_OBJECT<TV,T_MESH>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const // number_nodes must be set elsewhere
{
    return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Append_Particles_And_Create_Copy(*this,new_particles,particle_indices);
}
//#####################################################################
// Function Discard_Valence_Zero_Particles_And_Renumber
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Discard_Valence_Zero_Particles_And_Renumber(ARRAY<int>& condensation_mapping)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Discard_Valence_Zero_Particles_And_Renumber(*this,condensation_mapping);
}
//#####################################################################
// Function Union_Mesh_Objects_Relatively
//#####################################################################
template<class TV,class T_MESH> typename MESH_OBJECT<TV,T_MESH>::T_DERIVED_OBJECT* MESH_OBJECT<TV,T_MESH>::
Union_Mesh_Objects_Relatively(const ARRAY<T_DERIVED_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames)
{
    T_DERIVED_OBJECT* object=Create();
    Union_Mesh_Objects_Relatively(object,object_list,relative_frames);
    return object;
}
//#####################################################################
// Function Union_Mesh_Objects_Relatively
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Union_Mesh_Objects_Relatively(T_DERIVED_OBJECT *object,const ARRAY<T_DERIVED_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Union_Mesh_Objects_Relatively(object,object_list,relative_frames);
}
//#####################################################################
// Function Mark_Nodes_Referenced
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const
{
    mesh.Mark_Nodes_Referenced(marks,mark);
}
//#####################################################################
// Function Volume
//#####################################################################
struct UNUSABLE{};
template<class T_SIMPLICIAL_OBJECT> typename T_SIMPLICIAL_OBJECT::SCALAR
Filled_Volume_Helper(const T_SIMPLICIAL_OBJECT& object, typename DISABLE_IF<(T_SIMPLICIAL_OBJECT::MESH::dimension==(T_SIMPLICIAL_OBJECT::VECTOR_T::m-1)),UNUSABLE>::TYPE unused=UNUSABLE())
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class T_SIMPLICIAL_OBJECT> typename T_SIMPLICIAL_OBJECT::SCALAR
Filled_Volume_Helper(const T_SIMPLICIAL_OBJECT& object, typename ENABLE_IF<(T_SIMPLICIAL_OBJECT::MESH::dimension==(T_SIMPLICIAL_OBJECT::VECTOR_T::m-1)),UNUSABLE>::TYPE unused=UNUSABLE())
{
    typedef typename T_SIMPLICIAL_OBJECT::SCALAR T;typedef typename T_SIMPLICIAL_OBJECT::VECTOR_T TV;//enum WORKAROUND{d=T_SIMPLICIAL_OBJECT::MESH::dimension};
    static const int d=T_SIMPLICIAL_OBJECT::MESH::dimension;
    if(d!=TV::m-1) PHYSBAM_FATAL_ERROR("only codimension 1 objects can be filled");
    const TV base=object.particles.X(object.mesh.elements(1)[1]);
    T scaled_volume=0; // (d+1)!*volume
    for(int t=1;t<=object.mesh.elements.m;t++){const VECTOR<int,d+1>& nodes=object.mesh.elements(t);
        MATRIX<T,TV::m,d+1> DX;for(int i=1;i<=nodes.m;i++) DX.Column(i)=object.particles.X(nodes[i])-base;
        scaled_volume+=DX.Parallelepiped_Measure();}
    return (T)1/FACTORIAL<d+1>::value*scaled_volume;
}

template<class T> T
Filled_Volume_Helper(const MESH_OBJECT<VECTOR<T,1>, POINT_SIMPLEX_MESH>& object)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV,class T_MESH>
typename TV::SCALAR MESH_OBJECT<TV,T_MESH>::Volumetric_Volume()
{
    if(!mesh.elements.m) PHYSBAM_FATAL_ERROR("mesh has no elements");
    return Filled_Volume_Helper(*this);
}
//#####################################################################
template class MESH_OBJECT<VECTOR<float,1>,POINT_SIMPLEX_MESH>;
template class MESH_OBJECT<VECTOR<float,1>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<float,2>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<float,2>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,TETRAHEDRON_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,HEXAHEDRON_MESH>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MESH_OBJECT<VECTOR<double,1>,POINT_SIMPLEX_MESH>;
template class MESH_OBJECT<VECTOR<double,1>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<double,2>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<double,2>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,TETRAHEDRON_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,HEXAHEDRON_MESH>;
#endif
