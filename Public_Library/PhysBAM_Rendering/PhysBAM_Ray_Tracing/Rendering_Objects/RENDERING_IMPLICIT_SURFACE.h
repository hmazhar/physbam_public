//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_IMPLICIT_SURFACE
//#####################################################################
#ifndef __RENDERING_IMPLICIT_SURFACE__
#define __RENDERING_IMPLICIT_SURFACE__

#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T>
class RENDERING_IMPLICIT_SURFACE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::inverse_transform;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    IMPLICIT_OBJECT<TV>* implicit_surface;

    RENDERING_IMPLICIT_SURFACE(IMPLICIT_OBJECT<TV>* implicit_surface_input)
        :implicit_surface(implicit_surface_input)
    {}

    template<class T_GRID>
    RENDERING_IMPLICIT_SURFACE(T_GRID& grid_input,typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR& phi_input)
    {
        implicit_surface=new typename LEVELSET_POLICY<T_GRID>::LEVELSET_IMPLICIT_OBJECT(grid_input,phi_input);
    }

    virtual ~RENDERING_IMPLICIT_SURFACE()
    {delete implicit_surface;}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(implicit_surface->Intersection(object_space_ray,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(implicit_surface->Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return implicit_surface->Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return implicit_surface->Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return implicit_surface->Boundary(Object_Space_Point(location),small_number);}

    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {return (*implicit_surface)(Object_Space_Point(location));}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {TRIANGULATED_SURFACE<T>* surface=TESSELLATION::Generate_Triangles(*implicit_surface);surface->Update_Triangle_List();return surface;}

    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {return implicit_surface->box;}

//#####################################################################
};
}
#endif
