//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_SPHERE_PARTITION
//##################################################################### 
#ifndef __RENDERING_SPHERE_PARTITION__
#define __RENDERING_SPHERE_PARTITION__

#include <PhysBAM_Geometry/Spatial_Acceleration/SPHERE_PARTITION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_SPHERE_PARTITION:public RENDERING_OBJECT<T>
{
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::material;

    SPHERE_PARTITION<T> sphere_partition;
    
    RENDERING_SPHERE_PARTITION(const int number_input)
        :sphere_partition(number_input)
    {}

    virtual ~RENDERING_SPHERE_PARTITION()
    {}

    bool Intersection(RAY<VECTOR<T,3> >& ray) const PHYSBAM_OVERRIDE
    {RAY<VECTOR<T,3> > object_space_ray=Object_Space_Ray(ray);
    if(sphere_partition.Intersection(object_space_ray,small_number)){ray.semi_infinite=false;
        if(material.volumetric && sphere_partition.Inside(object_space_ray.Point((T).5*object_space_ray.t_max),small_number)){
            T volumetric_integration_step=material.Volumetric_Integration_Step();
            if(object_space_ray.t_max > volumetric_integration_step){ray.t_max=volumetric_integration_step;ray.aggregate_id=-1;return true;}}   
        ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    VECTOR<T,3> Normal(const VECTOR<T,3>& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(sphere_partition.Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const VECTOR<T,3>& location) const  PHYSBAM_OVERRIDE
    {return sphere_partition.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return sphere_partition.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return sphere_partition.Boundary(Object_Space_Point(location),small_number);}

//#####################################################################
};    
}
#endif

