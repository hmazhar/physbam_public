//#####################################################################
// Copyright 2006-2008, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY
//#####################################################################
#ifndef __RIGID_GEOMETRY__
#define __RIGID_GEOMETRY__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_STATE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>

namespace PhysBAM{

template<class TV>
class RIGID_GEOMETRY:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    typedef typename IF<TV::m==2,SEGMENT_HIERARCHY<TV>,TRIANGLE_HIERARCHY<T> >::TYPE T_SIMPLEX_HIERARCHY;

public:
    IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_GEOMETRY<TV> >* implicit_object; // implicit representation of geometry
    T_SIMPLICIAL_OBJECT* simplicial_object; // discrete representation of geometry
    RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection;
    int particle_index;
    bool is_static; // not saved to file - indicates whether this object is static in the scene
    T surface_roughness; // small number indicating errors in the geometry
    T coefficient_of_friction;

    T_ORIENTED_BOX oriented_box;
    RANGE<TV> axis_aligned_bounding_box;
    bool bounding_box_up_to_date;
    std::string name; // not saved to file - for convenience and debugging.
    T_SIMPLEX_HIERARCHY* moving_simplex_hierarchy; // bounds moving simplices in world space, to accelerate crossover check

    ARRAY<STRUCTURE<TV>*> structures;
    COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>* impulse_accumulator;

    RIGID_GEOMETRY(RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,bool create_collision_geometry,int index=0);
    virtual ~RIGID_GEOMETRY();

    template<class T_STRUCTURE>
    T_STRUCTURE Find_Structure(const int index=1) const
    {return Find_Type<T_STRUCTURE>(structures,index);}

    void Set_Name(const std::string& name_input)
    {name=name_input;}

    virtual std::string Name() const {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("RIGID_GEOMETRY<VECTOR<T,%d> >",TV::m);}

    void Set_Surface_Roughness(const T surface_roughness_input=(T)1e-6)
    {surface_roughness=surface_roughness_input;}

    void Set_Coefficient_Of_Friction(const T coefficient_input=.5)
    {coefficient_of_friction=coefficient_input;}

    static T Coefficient_Of_Friction(const RIGID_GEOMETRY<TV>& body1,const RIGID_GEOMETRY<TV>& body2)
    {return min(body1.coefficient_of_friction,body2.coefficient_of_friction);}

    RAY<TV> Object_Space_Ray(const RAY<TV>& world_space_ray) const
    {RAY<TV> transformed_ray(Object_Space_Point(world_space_ray.endpoint),Object_Space_Vector(world_space_ray.direction));
    transformed_ray.semi_infinite=world_space_ray.semi_infinite;transformed_ray.t_max=world_space_ray.t_max;
    transformed_ray.aggregate_id=world_space_ray.aggregate_id;
    return transformed_ray;}

    TV Object_Space_Point(const TV& world_space_point) const
    {return Frame().Inverse_Times(world_space_point);}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return Rotation().Inverse_Rotate(world_space_vector);}

    TV World_Space_Point(const TV& object_space_point) const
    {return Frame()*object_space_point;}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return Rotation().Rotate(object_space_vector);}

    TV& X() PHYSBAM_ALWAYS_INLINE
    {return rigid_geometry_collection.particles.X(particle_index);}
    
    const TV& X() const PHYSBAM_ALWAYS_INLINE
    {return rigid_geometry_collection.particles.X(particle_index);}

    ROTATION<TV>& Rotation() PHYSBAM_ALWAYS_INLINE
    {return rigid_geometry_collection.particles.rotation(particle_index);}
    
    const ROTATION<TV>& Rotation() const PHYSBAM_ALWAYS_INLINE
    {return rigid_geometry_collection.particles.rotation(particle_index);}

    const FRAME<TV> Frame() const PHYSBAM_ALWAYS_INLINE
    {return FRAME<TV>(X(),Rotation());}

    void Set_Frame(const FRAME<TV>& frame) PHYSBAM_ALWAYS_INLINE
    {X()=frame.t;Rotation()=frame.r;}

    TWIST<TV>& Twist() PHYSBAM_ALWAYS_INLINE
    {return rigid_geometry_collection.particles.twist(particle_index);}

    const TWIST<TV>& Twist() const PHYSBAM_ALWAYS_INLINE
    {return rigid_geometry_collection.particles.twist(particle_index);}

    static TV Pointwise_Object_Velocity(const TWIST<TV>& twist,const TV& t,const TV& X)
    {return twist.linear+TV::Cross_Product(twist.angular,X-t);}

    static TV Relative_Velocity(const RIGID_GEOMETRY<TV>& geometry1,const RIGID_GEOMETRY<TV>& geometry2,const TV& world_point)
    {return geometry1.Pointwise_Object_Velocity(world_point)-geometry2.Pointwise_Object_Velocity(world_point);}

    static TV Relative_Velocity_At_Geometry1_Particle(const RIGID_GEOMETRY<TV>& geometry1,const RIGID_GEOMETRY<TV>& geometry2,const TV& world_point,const int particle_index)
    {return geometry1.Pointwise_Object_Velocity_At_Particle(world_point,particle_index)-geometry2.Pointwise_Object_Velocity(world_point);}

    static T_SPIN Relative_Angular_Velocity(const RIGID_GEOMETRY<TV>& geometry1,const RIGID_GEOMETRY<TV>& geometry2) // make sure the angular velocities are updated before calling this!
    {return geometry1.Twist().angular-geometry2.Twist().angular;}

    static TWIST<TV> Relative_Twist(const RIGID_GEOMETRY<TV>& geometry1,const RIGID_GEOMETRY<TV>& geometry2,const TV& world_point)
    {return TWIST<TV>(geometry1.Pointwise_Object_Velocity(world_point)-geometry2.Pointwise_Object_Velocity(world_point),geometry1.Twist().angular-geometry2.Twist().angular);}

    const T_ORIENTED_BOX& Oriented_Bounding_Box() const
    {assert(bounding_box_up_to_date);return oriented_box;}

    const RANGE<TV>& Axis_Aligned_Bounding_Box() const
    {assert(bounding_box_up_to_date);return axis_aligned_bounding_box;}

    bool Simplex_Inside(const TV& location,const T collision_thickness) const
    {return simplicial_object->Inside(Object_Space_Point(location),collision_thickness);}

    bool Simplex_Outside(const TV& location,const T collision_thickness) const
    {return simplicial_object->Outside(Object_Space_Point(location),collision_thickness);}

    bool Bounding_Boxes_Intersect(const RIGID_GEOMETRY<TV>& rigid_geometry,const T thickness=0) const
    {if(!thickness) return Axis_Aligned_Bounding_Box().Lazy_Intersection(rigid_geometry.Axis_Aligned_Bounding_Box()) // check axis aligned first for speed
        && Oriented_Bounding_Box().Intersection(rigid_geometry.Oriented_Bounding_Box());
    return Axis_Aligned_Bounding_Box().Intersection(rigid_geometry.Axis_Aligned_Bounding_Box(),thickness)
        && Oriented_Bounding_Box().Thickened(thickness).Intersection(rigid_geometry.Oriented_Bounding_Box());} // Thickened is rather expensive; avoid it where possible

    TV Simplex_Normal(const TV& location,const int aggregate=0) const
    {return World_Space_Vector(simplicial_object->Normal(Object_Space_Point(location),aggregate));}

    bool Simplex_Boundary(const TV& location,const T collision_thickness) const
    {return simplicial_object->Boundary(Object_Space_Point(location),collision_thickness);}    

//#####################################################################
    void Update_Bounding_Box();
    void Update_Bounding_Box_From_Implicit_Geometry();
    void Add_Structure(STRUCTURE<TV>& structure); // set up acceleration stuctures for certain types of structures
    static void Print_Names(const RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection,const ARRAY<int>& indices);
    virtual TV Pointwise_Object_Velocity(const TV& X) const;
    virtual TV Pointwise_Object_Velocity_At_Particle(const TV& X,const int particle_index) const;
    TV Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const;
    TV Implicit_Geometry_Normal(const TV& location,const int aggregate=-1) const;
    TV Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const;
    bool Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value=0) const;
    bool Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value=0) const;
    T Implicit_Geometry_Extended_Value(const TV& location) const;
    const RANGE<TV>& Object_Space_Bounding_Box() const;
    bool Simplex_Intersection(RAY<TV>& ray,const T collision_thickness) const;
    RANGE<TV> World_Space_Simplex_Bounding_Box(const int id) const;
    typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX World_Space_Simplex(const int id) const;
    RANGE<TV> World_Space_Simplex_Bounding_Box(const int id,const FRAME<TV>& frame) const;
    void Interpolate_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,const T time,RIGID_GEOMETRY_STATE<TV>& interpolated_state);
    void Compute_Velocity_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,RIGID_GEOMETRY_STATE<TV>& result_state);
protected:
    void Remove_Structure(STRUCTURE<TV>* structure);
    template<class T> friend void Add_Structure_Helper(RIGID_GEOMETRY<VECTOR<T,1> >& self,STRUCTURE<VECTOR<T,1> >& structure);
    template<class T> friend void Add_Structure_Helper(RIGID_GEOMETRY<VECTOR<T,2> >& self,STRUCTURE<VECTOR<T,2> >& structure);
    template<class T> friend void Add_Structure_Helper(RIGID_GEOMETRY<VECTOR<T,3> >& self,STRUCTURE<VECTOR<T,3> >& structure);
};
template<class TV,class SCALAR> struct REPLACE_FLOATING_POINT<RIGID_GEOMETRY<TV>,SCALAR>{typedef RIGID_GEOMETRY<typename REPLACE_FLOATING_POINT<TV,SCALAR>::TYPE> TYPE;};
}
#endif
