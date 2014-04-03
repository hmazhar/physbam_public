#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_OCTREE_IMPLICIT_SURFACE
//##################################################################### 
#ifndef __RENDERING_OCTREE_IMPLICIT_SURFACE__
#define __RENDERING_OCTREE_IMPLICIT_SURFACE__

#include <PhysBAM_Geometry/Grids_Dyadic_Computations/DUALCONTOUR_OCTREE.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_OCTREE_IMPLICIT_SURFACE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    DYADIC_IMPLICIT_OBJECT<TV>* implicit_surface;

    RENDERING_OCTREE_IMPLICIT_SURFACE(const OCTREE_GRID<T>& octree_grid_input,const ARRAY<T>& phi_input,const int maximum_depth=0)
    {
        implicit_surface=new DYADIC_IMPLICIT_OBJECT<TV>(octree_grid_input,phi_input,maximum_depth);
    }
    
    virtual ~RENDERING_OCTREE_IMPLICIT_SURFACE()
    {delete implicit_surface;}

    LEVELSET_OCTREE<T>& Octree_Levelset()
    {return implicit_surface->octree_levelset;}
    
    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    bool intersection=implicit_surface->Intersection(object_space_ray,small_number);
    ray.semi_infinite=object_space_ray.semi_infinite;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return intersection;}

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

    TRIANGULATED_SURFACE<T>* Generate_Triangles()const PHYSBAM_OVERRIDE
    {DUALCONTOUR_OCTREE<T> contour(&implicit_surface->levelset);return contour.Get_Triangulated_Surface();}

//#####################################################################
};    
}
#endif
#endif
