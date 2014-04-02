//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_GEOMETRY_COLLECTION
//#####################################################################
#ifndef __COLLISION_GEOMETRY_COLLECTION__
#define __COLLISION_GEOMETRY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
namespace PhysBAM{

template<class TV> class RAY;

template<class TV>
class COLLISION_GEOMETRY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>* spatial_partition;
private:
    ARRAY<COLLISION_GEOMETRY_ID> deletion_list;
    ARRAY<bool,COLLISION_GEOMETRY_ID> owns_collision_geometry;
public:
    ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID> bodies;
    HASHTABLE<int,COLLISION_GEOMETRY_ID> geometry_id_to_collision_geometry_id;
    HASHTABLE<COLLISION_GEOMETRY_ID,int> collision_geometry_id_to_geometry_id;
    T collision_body_thickness;

    COLLISION_GEOMETRY_COLLECTION();
    virtual ~COLLISION_GEOMETRY_COLLECTION();

    COLLISION_GEOMETRY<TV>& operator()(const COLLISION_GEOMETRY_ID id)
    {return *bodies(id);}

    const COLLISION_GEOMETRY<TV>& operator()(const COLLISION_GEOMETRY_ID id) const
    {return *bodies(id);}

    COLLISION_GEOMETRY_ID Size() const
    {return bodies.Size();}

    bool Is_Active(COLLISION_GEOMETRY_ID id) const
    {return (bodies(id)!=0 && bodies(id)->active);}

    void Set_Collision_Body_Thickness(T collision_body_thickness_input)
    {
        collision_body_thickness=collision_body_thickness_input;if(spatial_partition) spatial_partition->Set_Collision_Body_Thickness(collision_body_thickness);
        for(COLLISION_GEOMETRY_ID id(1);id<=bodies.Size();++id) bodies(id)->Set_Collision_Thickness(collision_body_thickness);
    }

    COLLISION_GEOMETRY<TV>* Get_Collision_Geometry(int geometry_id)
    {return bodies(geometry_id_to_collision_geometry_id.Get(geometry_id));}
    
    bool Has_Collision_Body(int geometry_id)
    {return geometry_id_to_collision_geometry_id.Contains(geometry_id);}

//#####################################################################
    COLLISION_GEOMETRY_ID Add_Body(COLLISION_GEOMETRY<TV>* collision_body,const int geometry_id,bool owns_collision_geometry_input);
    void Add_Bodies(COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list);
    void Remove_Body(COLLISION_GEOMETRY_ID id);
    bool Intersection_Between_Points(const TV& x1,const TV& x2,COLLISION_GEOMETRY_ID& body_id,int& triangle_id,TV& intersection_point,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Intersection_Between_Points(const TV& x1,const TV& x2,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Intersection_With_Any_Simplicial_Object(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Closest_Non_Intersecting_Point_Of_Any_Simplicial_Object(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Inside_Any_Simplex_Of_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Inside_Any_Body(const TV& location,const T thickness_over_two,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Implicit_Geometry_Lazy_Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    TV Closest_Boundary_Point(const TV& location,const T max_distance,T& distance,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    void Update_Spatial_Partition(const SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC heuristic,const int number_of_boxes,const T voxel_size_scale_factor);
    void Update_Bounding_Boxes();
//#####################################################################
};
}
#endif
