//#####################################################################
// Copyright 2006-2009, Nipun Kwatra, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_COLLECTION
//#####################################################################
#ifndef __RIGID_GEOMETRY_COLLECTION__
#define __RIGID_GEOMETRY_COLLECTION__

#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>

namespace PhysBAM{

template<class TV> struct ALLOCATE_HELPER{virtual RIGID_GEOMETRY<TV>* Create(int index=0)=0;virtual ~ALLOCATE_HELPER(){}};
template<class TV,class ID> class STRUCTURE_LIST;
template<class TV> class STRUCTURE;

template<class TV>
class RIGID_GEOMETRY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    RIGID_GEOMETRY_PARTICLES<TV>& particles;
    COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list;
    STRUCTURE_LIST<TV,int>& structure_list;
    bool always_create_structure;
    HASHTABLE<std::string,int>& structure_hash; // maps to id
    mutable bool is_stale_key,is_stale_active;
    mutable ARRAY<int> *frame_list_key,*frame_list_active;
    mutable ARRAY<std::string,int> rigid_body_names;
    bool check_stale;
    int last_read_key,last_read_active;
    ALLOCATE_HELPER<TV>* allocate_helper;
    RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>* rigid_geometry_example_velocities;    
    ARRAY<int> static_rigid_geometry,kinematic_rigid_geometry;
    bool owns_particles;

    RIGID_GEOMETRY_COLLECTION(RIGID_GEOMETRY_PARTICLES<TV>& particles_input,RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>* rigid_geometry_example_velocities_input,
        COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input=0,ALLOCATE_HELPER<TV>* allocate_helper_input=0);
    RIGID_GEOMETRY_COLLECTION(RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>* rigid_geometry_example_velocities_input,COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input=0,ALLOCATE_HELPER<TV>* allocate_helper_input=0);
    virtual ~RIGID_GEOMETRY_COLLECTION();
    
    RIGID_GEOMETRY<TV>& Rigid_Geometry(const int particle_index)
    {return *particles.rigid_geometry(particle_index);}

    const RIGID_GEOMETRY<TV>& Rigid_Geometry(const int particle_index) const
    {return *particles.rigid_geometry(particle_index);}

    void Deactivate_Geometry(const int p)
    {assert(Exists(p) && Rigid_Geometry(p).particle_index>0);Rigid_Geometry(p).particle_index=-Rigid_Geometry(p).particle_index;}

    void Reactivate_Geometry(const int p)
    {assert(Exists(p) && Rigid_Geometry(p).particle_index<0);Rigid_Geometry(p).particle_index=-Rigid_Geometry(p).particle_index;}

    RIGID_GEOMETRY<TV>* New_Body(int index);

//#####################################################################
    bool Exists(const int particle) const;
    bool Is_Active(const int particle) const;
    void Update_Kinematic_Particles();
    int Add_Rigid_Geometry(STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file);
    int Add_Rigid_Geometry(RIGID_GEOMETRY<TV>* rigid_geometry,STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,
        const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file);
    bool Register_Analytic_Replacement_Structure(const std::string& filename,const T scaling_factor,STRUCTURE<TV>* structure); // passing in zero skips reading
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    bool Find_Or_Read_Structure(const STREAM_TYPE stream_type,ARRAY<int>& structure_id,const std::string& filename,const T scaling_factor,const TV& center);
#endif
    void Destroy_Unreferenced_Geometry();
//#####################################################################
};
}
#endif
