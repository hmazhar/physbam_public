//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_VOXELS
//#####################################################################
#ifndef __RENDERING_VOXELS__
#define __RENDERING_VOXELS__

#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_VOXELS:public RENDERING_OBJECT<T>
{
    using RENDERING_OBJECT<T>::Inside;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
    
    BOX<VECTOR<T,3> > box; // box containing the voxelized data
    bool precompute_single_scattering;

    RENDERING_VOXELS()
        :precompute_single_scattering(false)
    {}

    virtual ~RENDERING_VOXELS()
    {}

    bool Intersection(RAY<VECTOR<T,3> >& ray)  const PHYSBAM_OVERRIDE
    {RAY<VECTOR<T,3> > object_space_ray=Object_Space_Ray(ray);
    bool box_intersection=false;
    if(INTERSECTION::Intersects(object_space_ray,box)){box_intersection=true;ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;}
    return box_intersection;}

    bool Get_Intersection_Range(const RAY<VECTOR<T,3> >& ray,T& start_t,T& end_t) const PHYSBAM_OVERRIDE
    {return INTERSECTION::Get_Intersection_Range(Object_Space_Ray(ray),box,start_t,end_t);}

    VECTOR<T,3> Normal(const VECTOR<T,3>& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(box.Normal(aggregate));}

    bool Inside(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return box.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return box.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return box.Boundary(Object_Space_Point(location),small_number);}

    void Preprocess_Efficiency_Structures(RENDER_WORLD<T>& world) PHYSBAM_OVERRIDE
    {if(precompute_single_scattering) Precompute_Light_Data(true,world);}

//#####################################################################
    void Precompute_Light_Data(bool use_fast_precomputation,RENDER_WORLD<T>& world);
    virtual T Volumetric_Integration_Step(const RAY<VECTOR<T,3> >& ray, const T xi) const=0;
    virtual bool Use_Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const {return false;}
    virtual void Set_Precomputed_Light_Data(const int location_index,const int light_index,const VECTOR<T,3>& light_value){}
    virtual void Set_Precomputed_Light_Valid(const int location_index,const int light_index,const bool value){}
    virtual VECTOR<T,3> Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const {return VECTOR<T,3>(0,0,0);}
    virtual T Source_Term(const int source_term_index,const VECTOR<T,3>& location) const {return 0;}
    virtual void Get_Node_Locations(ARRAY<VECTOR<T,3> >& locations){}
protected:
    virtual void Postprocess_Light_Field(){}
    virtual void Prepare_For_Precomputation(RENDER_WORLD<T>& world){};
//#####################################################################
};   
}
#endif

