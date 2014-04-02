//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Jon Gretarsson, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAY
//#####################################################################
#ifndef __RAY__
#define __RAY__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <cmath>
namespace PhysBAM{

template<class TV>
class RAY
{
public:
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    TV endpoint; // endpoint of the ray where t=0
    TV direction; // direction the ray sweeps out - unit vector
    bool semi_infinite; // indicates whether the ray is semi_infinite or should stop at t_max
    T t_max; // maximum value of t allowed for the ray
    int aggregate_id; // indicates the piece of an aggregate object that is intersected by t_max
    enum LOCATION {START_POINT,END_POINT,INTERIOR_POINT,LOCATION_UNKNOWN};
    LOCATION intersection_location; // indicates what type of intersection happened, LOCATION_UNKNOWN is used if not computed
    RANGE<TV> bounding_box;

private:
    // used for triangle hierarchy fast lazy_box_intersection 
    bool computed_lazy_box_intersection_acceleration_data;
public:
    TV inverse_direction;
    TV_INT direction_is_negative;

    RAY()
        :endpoint(TV()),semi_infinite(true),t_max(0),aggregate_id(0),intersection_location(LOCATION_UNKNOWN),computed_lazy_box_intersection_acceleration_data(false)
    {
         direction=TV::Axis_Vector(TV::dimension);
    }

    RAY(const TV& endpoint_input,const TV& direction_input,const bool already_normalized=false)
        :endpoint(endpoint_input),direction(direction_input),semi_infinite(true),t_max(0),aggregate_id(0),intersection_location(LOCATION_UNKNOWN),computed_lazy_box_intersection_acceleration_data(false)
    {if(!already_normalized) direction.Normalize();}

    RAY(const SEGMENT_1D<T>& segment)
        :endpoint(segment.x1),direction(segment.x2-segment.x1),semi_infinite(false),aggregate_id(0),intersection_location(LOCATION_UNKNOWN),computed_lazy_box_intersection_acceleration_data(false)
    {
        STATIC_ASSERT(TV::dimension==1);
        t_max=direction.Normalize();
    }

    RAY(const SEGMENT_2D<T>& segment)
        :endpoint(segment.x1),direction(segment.x2-segment.x1),semi_infinite(false),aggregate_id(0),intersection_location(LOCATION_UNKNOWN),computed_lazy_box_intersection_acceleration_data(false)
    {
        STATIC_ASSERT(TV::dimension==2);
        t_max=direction.Normalize();
    }

    RAY(const SEGMENT_3D<T>& segment)
        :endpoint(segment.x1),direction(segment.x2-segment.x1),semi_infinite(false),aggregate_id(0),intersection_location(LOCATION_UNKNOWN),computed_lazy_box_intersection_acceleration_data(false)
    {
        STATIC_ASSERT(TV::dimension==3);
        t_max=direction.Normalize();
    }

    void Initialize(const TV& endpoint_input,const TV& direction_input,const bool already_normalized=false)
    {endpoint=endpoint_input;direction=direction_input;semi_infinite=true;t_max=0;aggregate_id=0;intersection_location=LOCATION_UNKNOWN;computed_lazy_box_intersection_acceleration_data=false;
    if(!already_normalized) direction.Normalize();}

    void Save_Intersection_Information(RAY<TV>& storage_ray) const
    {storage_ray.semi_infinite=semi_infinite;storage_ray.t_max=t_max;storage_ray.aggregate_id=aggregate_id;storage_ray.intersection_location=intersection_location;}

    void Restore_Intersection_Information(const RAY<TV>& storage_ray)
    {semi_infinite=storage_ray.semi_infinite;t_max=storage_ray.t_max;aggregate_id=storage_ray.aggregate_id;intersection_location=storage_ray.intersection_location;}

    TV Point(const T t) const // finds the point on the ray, given by the parameter t
    {return endpoint+t*direction;}

    void Compute_Bounding_Box()
    {assert(!semi_infinite);bounding_box=RANGE<TV>::Bounding_Box(endpoint,Point(t_max));}

    T Parameter(const TV& point) const // finds the parameter t, given a point that lies on the ray
    {int axis=direction.Dominant_Axis();return (point[axis]-endpoint[axis])/direction[axis];}

    TV Reflected_Direction(const TV& normal) const
    {return 2*TV::Dot_Product(-direction,normal)*normal+direction;}

    static bool Create_Non_Degenerate_Ray(const TV& endpoint,const TV& length_and_direction,RAY<TV>& ray)
    {T length_squared=length_and_direction.Magnitude_Squared();
    if(length_squared>0){ray.t_max=sqrt(length_squared);ray.endpoint=endpoint;ray.direction=length_and_direction/ray.t_max;ray.semi_infinite=false;return true;}
    else return false;}

//#####################################################################
    void Compute_Lazy_Box_Intersection_Acceleration_Data();
//#####################################################################
};

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class TV> std::ostream &operator<<(std::ostream &output,const RAY<TV> &ray)
{output<<"endpoint = "<<ray.endpoint<<", direction = "<<ray.direction<<", ";
if(ray.semi_infinite) output<<"semi infinite";else output<<"t_max = "<<ray.t_max;output<<std::endl;
return output;
}
#endif

}
#endif
