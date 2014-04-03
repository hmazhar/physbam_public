//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_SHOCKS
//#####################################################################
#ifndef __RENDERING_SHOCKS__
#define __RENDERING_SHOCKS__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T>
class RENDERING_SHOCKS:public RENDERING_VOXELS<T>
{
    typedef VECTOR<T,3> TV;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
public:
    using RENDERING_VOXELS<T>::box;using RENDERING_VOXELS<T>::Intersection; // silence -Woverloaded-virtual

    const GRID<TV>& grid;
    const T_ARRAYS_SCALAR& density;
    const T_ARRAYS_SCALAR& pressure;
    T gradient_threshold;
    T refraction_multiplier;
    T volumetric_step,fine_volumetric_step;
    T skip_next_intersection_factor;
    T dale_constant;
    bool use_pressure_for_intersection,use_pressure_for_rarefaction;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
    mutable TV last_intersection_point;
    mutable T last_intersection_t;

    RENDERING_SHOCKS(const GRID<TV>& grid_input,const T_ARRAYS_SCALAR& density_input,const T_ARRAYS_SCALAR& pressure_input,
        const T gradient_threshold_input,const T refraction_multiplier_input,const T volumetric_step_input,
        const T fine_volumetric_step_input,const T skip_next_intersection_factor_input,
        const bool use_pressure_for_intersection_input,const bool use_pressure_for_rarefaction_input)
        :grid(grid_input),density(density_input),pressure(pressure_input),gradient_threshold(gradient_threshold_input),
        refraction_multiplier(refraction_multiplier_input),volumetric_step(volumetric_step_input),
        fine_volumetric_step(fine_volumetric_step_input),skip_next_intersection_factor(skip_next_intersection_factor_input),
        use_pressure_for_intersection(use_pressure_for_intersection_input),
        use_pressure_for_rarefaction(use_pressure_for_rarefaction_input)
    {box=grid.Domain();
    dale_constant=2.26e-4;
    LOG::cout<<"RENDERING_SHOCKS"<<std::endl;
    LOG::cout<<"gradient_threshold="<<gradient_threshold<<", refraction_multiplier="<<refraction_multiplier<<std::endl;
    LOG::cout<<"volumetric_step="<<volumetric_step<<", fine_volumetric_step="<<fine_volumetric_step<<", skip_next_intersection_factor="<<skip_next_intersection_factor<<std::endl;
    LOG::cout<<"use_pressure_for_intersection="<<use_pressure_for_intersection<<", use_pressure_for_rarefaction="<<use_pressure_for_rarefaction<<std::endl;
    }

    virtual ~RENDERING_SHOCKS()
    {}

    virtual T Volumetric_Integration_Step(const RAY<VECTOR<T,3> > &ray,const T xi) const PHYSBAM_OVERRIDE
    {return xi*volumetric_step;}

    virtual T Index_Of_Refraction(const TV& world_space_location) const PHYSBAM_OVERRIDE
    {TV object_space_location=Object_Space_Point(world_space_location);
    return (T)1+dale_constant*refraction_multiplier*Medium_For_Rarefaction(object_space_location);}

    T Density(const TV& location) const
    {return interpolation.Clamped_To_Array(grid,density,location);}

    T Pressure(const TV& location) const
    {return interpolation.Clamped_To_Array(grid,pressure,location);}

    T Medium_For_Intersection(const TV& location) const
    {if(use_pressure_for_intersection) return Pressure(location);
    else return Density(location);}

    T Medium_For_Rarefaction(const TV& location) const
    {if(use_pressure_for_rarefaction) return Pressure(location);
    else return Density(location);}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray); // TODO: do we need this?
    T start_t,end_t;
    if(!INTERSECTION::Get_Intersection_Range(object_space_ray,box,start_t,end_t)) return false;

    if(start_t==0) start_t+=volumetric_step*skip_next_intersection_factor;

    T current_t=start_t;
    bool last_segment=false;
    int step_count=0;
    while(!last_segment){
        step_count++;
        TV current_start_position=object_space_ray.Point(current_t);
        T current_volumetric_step=volumetric_step;
        if(current_t+current_volumetric_step>end_t){last_segment=true;current_volumetric_step=end_t-current_t;}
        TV current_end_position=object_space_ray.Point(current_t+current_volumetric_step);

        T medium_start=Medium_For_Intersection(current_start_position);
        T medium_end=Medium_For_Intersection(current_end_position);

        T current_gradient=abs(medium_end-medium_start)/current_volumetric_step;
        if(current_gradient>gradient_threshold){
            // refine to find more precise location of highest medium gradient.

            T max_gradient=current_gradient;
            T max_gradient_t=current_t+current_volumetric_step*(T).5;
            T end_max_grad_t=min(end_t,current_t+2*current_volumetric_step);
            bool last_max_gradient_segment=fine_volumetric_step>=current_volumetric_step; // skip if fine step larger than current
            while(!last_max_gradient_segment){
                T current_max_grad_step=fine_volumetric_step;
                if(current_t+current_max_grad_step>end_max_grad_t){
                    last_max_gradient_segment=true;current_max_grad_step=end_max_grad_t-current_t;}

                current_start_position=object_space_ray.Point(current_t);
                current_end_position=object_space_ray.Point(current_t+current_max_grad_step);
                medium_start=Medium_For_Intersection(current_start_position);
                medium_end=Medium_For_Intersection(current_end_position);

                current_gradient=abs(medium_end-medium_start)/current_max_grad_step;
                if(current_gradient > max_gradient){
                    max_gradient=current_gradient;
                    max_gradient_t=current_t+current_max_grad_step*(T).5;}
                current_t+=current_max_grad_step;}

            last_intersection_point=object_space_ray.Point(max_gradient_t);
            last_intersection_t=max_gradient_t;
            ray.semi_infinite=false;ray.t_max=max_gradient_t;return true;}
        current_t+=current_volumetric_step;}

    return false;}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {TV object_space_location=Object_Space_Point(location);

    TV dX=grid.DX();
    TV object_space_normal=TV(
    (Medium_For_Intersection(object_space_location+dX(1)*TV::Axis_Vector(1))-Medium_For_Intersection(object_space_location-dX(1)*TV::Axis_Vector(1)))/(2*dX(1)),
    (Medium_For_Intersection(object_space_location+dX(2)*TV::Axis_Vector(2))-Medium_For_Intersection(object_space_location-dX(2)*TV::Axis_Vector(2)))/(2*dX(2)),
    (Medium_For_Intersection(object_space_location+dX(3)*TV::Axis_Vector(3))-Medium_For_Intersection(object_space_location-dX(3)*TV::Axis_Vector(3)))/(2*dX(3))).Normalized();

    return World_Space_Vector(object_space_normal);}

    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {return box;}
//#####################################################################
};
}
#endif
