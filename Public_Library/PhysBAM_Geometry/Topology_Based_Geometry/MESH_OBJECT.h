//#####################################################################
// Copyright 2002-2009, Ron Fedkiw, Eilene Hao, Geoffrey Irving, Michael Lentine, Neil Molino, Joseph Teran, Robert Bridson, Andrew Selle, Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MESH_OBJECT
//##################################################################### 
#ifndef __MESH_OBJECT__
#define __MESH_OBJECT__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
namespace PhysBAM{
template<class TV,class T_MESH> struct MESH_TO_OBJECT; // mapping from mesh to most derived child of MESH_OBJECT
template<class TV> struct MESH_TO_OBJECT<TV,POINT_SIMPLEX_MESH>{typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::POINT_SIMPLICES TYPE;};
template<class TV> struct MESH_TO_OBJECT<TV,SEGMENT_MESH>{typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE TYPE;};
template<class TV> struct MESH_TO_OBJECT<TV,TRIANGLE_MESH>{typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT TYPE;};
template<class T> struct MESH_TO_OBJECT<VECTOR<T,3>,TETRAHEDRON_MESH>{typedef TETRAHEDRALIZED_VOLUME<T> TYPE;};
template<class T> struct MESH_TO_OBJECT<VECTOR<T,3>,HEXAHEDRON_MESH>{typedef HEXAHEDRALIZED_VOLUME<T> TYPE;};

template<class TV,class T_MESH>
class MESH_OBJECT:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    enum WORKAROUND1 {dimension=T_MESH::dimension};
    typedef T_MESH MESH;
    typedef typename MESH_TO_OBJECT<TV,T_MESH>::TYPE T_DERIVED_OBJECT;

    T_MESH& mesh;
    GEOMETRY_PARTICLES<TV>& particles;
    BOX<TV>* bounding_box;
    PARTICLE_PARTITION<TV>* particle_partition;
    TV_INT desired_particle_partition_counts;
private:
    bool need_destroy_mesh,need_destroy_particles;
protected:
    MESH_OBJECT(T_MESH& mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);
    virtual ~MESH_OBJECT();
public:

    void Discard_Valence_Zero_Particles_And_Renumber()
    {ARRAY<int> condensation_mapping;Discard_Valence_Zero_Particles_And_Renumber(condensation_mapping);}

    void Set_Desired_Particle_Partition_Size(const TV_INT& counts)
    {desired_particle_partition_counts=counts;}

//#####################################################################
    virtual void Clean_Memory();
    static T_DERIVED_OBJECT* Create();
    static T_DERIVED_OBJECT* Create(GEOMETRY_PARTICLES<TV>& particles);
    virtual void Refresh_Auxiliary_Structures() PHYSBAM_SEALED;
    void Update_Number_Nodes() PHYSBAM_OVERRIDE;
    void Update_Bounding_Box();
    void Initialize_Particle_Partition(const TV_INT& counts=TV_INT());
    STRUCTURE<TV>* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices=0) const PHYSBAM_OVERRIDE; // number_nodes must be set elsewhere
    void Discard_Valence_Zero_Particles_And_Renumber(ARRAY<int>& condensation_mapping);
    static T_DERIVED_OBJECT* Union_Mesh_Objects_Relatively(const ARRAY<T_DERIVED_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames);
    static void Union_Mesh_Objects_Relatively(T_DERIVED_OBJECT* object,const ARRAY<T_DERIVED_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames);
    void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const PHYSBAM_OVERRIDE;
    T Volumetric_Volume();
private:
    virtual void Refresh_Auxiliary_Structures_Helper()=0;
//#####################################################################
};
}
#endif
