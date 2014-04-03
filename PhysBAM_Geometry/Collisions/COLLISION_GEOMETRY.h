//#####################################################################
// Copyright 2005-2008, Elliot English, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_GEOMETRY
//#####################################################################
#ifndef __COLLISION_GEOMETRY__
#define __COLLISION_GEOMETRY__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/TYPED_STREAM.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> class BINDING_LIST;
template<class TV> class SOFT_BINDINGS;
template<class TV> class FRAME;
template<class TV> class RANGE;
template<class TV> class COLLISION_PARTICLE_STATE;

template<class TV>
class COLLISION_GEOMETRY:public NONCOPYABLE
{
private:
    typedef typename TV::SCALAR T;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_SIMPLEX;
    typedef typename IF<TV::dimension==2,T,typename IF<TV::dimension==1,ONE,TV>::TYPE>::TYPE T_WEIGHTS;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
public:
    typedef TV VECTOR_T;
    enum STATE_INDEX {FLUID_COLLISION_GEOMETRY_NEW_STATE=1,FLUID_COLLISION_GEOMETRY_OLD_STATE,FLUID_COLLISION_GEOMETRY_SAVED_OLD_STATE,FLUID_COLLISION_GEOMETRY_SAVED_NEW_STATE};

    bool active; // when false, body is ignored TODO: remove this
    T collision_thickness;
    bool add_to_spatial_partition;
    COLLISION_GEOMETRY_ID collision_geometry_id;
    COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>* impulse_accumulator;
    bool refine_nearby_fluid;

    ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* collision_geometries_for_rasterization;

    COLLISION_GEOMETRY();
    virtual ~COLLISION_GEOMETRY();

    void Set_Collision_Thickness(const T thickness=(T)1e-3)
    {collision_thickness=thickness;}

    void Set_Collision_Geometry_Id_Number(const COLLISION_GEOMETRY_ID collision_geometry_id_input)
    {collision_geometry_id=collision_geometry_id_input;}

//#####################################################################
    virtual TV Pointwise_Object_Velocity(const TV& X) const;
    virtual TV Pointwise_Object_Velocity(const int aggregate_id,const TV& X) const;
    virtual TV Pointwise_Object_Pseudo_Velocity(const int aggregate_id,const TV& X,const int state1,const int state2) const;
    // TODO: change name to Refresh_Auxiliary_Structures
    virtual void Update_Bounding_Box();
    virtual const RANGE<TV>& Axis_Aligned_Bounding_Box() const;
    virtual bool Has_Volumetric_Geometry() const;

    // interface independent
    virtual TV Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2,const int max_iterations,int* simplex_id,T* distance=0) const;

    // levelset interface
    virtual T Implicit_Geometry_Extended_Value(const TV& location) const;
    virtual TV Implicit_Geometry_Normal(const TV& location,const int aggregate=-1) const;
    // TODO: remove need for location_particle_index parameter
    virtual TV Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const;
    virtual TV Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const;
    virtual TV Implicit_Geometry_Closest_Point_On_Boundary(const TV& location,const T tolerance=0,const int max_iterations=1,T* distance=0) const;
    virtual bool Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value=0) const;
    virtual bool Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value=0) const;
    virtual bool Implicit_Geometry_Lazy_Inside_Extended_Levelset(const TV& location,T contour_value=0) const;

    // simplex interface
    virtual TV Simplex_World_Space_Point_From_Barycentric_Coordinates(const int simplex_id,const T_WEIGHTS& weights) const;
    virtual int Number_Of_Simplices() const;
    virtual bool Simplex_Intersection(RAY<TV>& ray) const;
    virtual bool Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const;
    virtual bool Inside_Any_Simplex(const TV& location,int& simplex_id) const;
    virtual bool Inside(const TV& location,const T thickness_over_two) const;
    virtual TV Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2=0,int* simplex_id=0,T* distance=0) const;
    virtual T_SIMPLEX World_Space_Simplex(const int simplex_id,const bool use_saved_state=false) const;
    virtual bool Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id) const;
    virtual bool Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id,
        POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const;
    virtual bool Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const;
    virtual void Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const;

    // push out
    virtual bool Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,T& hit_time,int& simplex_id,T& start_phi,T& end_phi,TV& end_body_normal,
         TV& body_velocity) const;
    virtual bool Push_Out_Point(TV& X,const T collision_distance,T& distance) const;

    // embedded collision
    virtual void Update_Intersection_Acceleration_Structures(const bool use_swept_triangle_hierarchy,const int state1=0,const int state2=0);

    virtual void Save_State(const int state_index,const T time=0);
    virtual void Restore_State(const int state_index);
    virtual void Average_States(const int state1, const int state2, const int result_state,const T interpolation_distance);
    virtual void Delete_State(const int state_index);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    virtual void Read_State(TYPED_ISTREAM& input,const int state_index);
    virtual void Write_State(TYPED_OSTREAM& output,const int state_index) const;
#endif
//#####################################################################
};
}
#endif
