#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_RLE_IMPLICIT_SURFACE
//#####################################################################
#ifndef __RENDERING_RLE_IMPLICIT_SURFACE__
#define __RENDERING_RLE_IMPLICIT_SURFACE__

#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_RLE_IMPLICIT_SURFACE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;

    RLE_IMPLICIT_OBJECT<TV> implicit_surface;

    RENDERING_RLE_IMPLICIT_SURFACE(RLE_GRID_3D<T>& grid_input,ARRAY<T>& phi_input)
        :implicit_surface(grid_input,phi_input)
    {
        PHYSBAM_FATAL_ERROR(); // this class is unnecessary now, but may be important if we need to speed up rle rendering
    }

    virtual ~RENDERING_RLE_IMPLICIT_SURFACE()
    {}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(implicit_surface.Intersection(object_space_ray,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    VECTOR<T,3> Normal(const VECTOR<T,3>& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(implicit_surface.Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return implicit_surface.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return implicit_surface.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return implicit_surface.Boundary(Object_Space_Point(location),small_number);}

    T Signed_Distance(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return implicit_surface(Object_Space_Point(location));}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {TRIANGULATED_SURFACE<T>* surface=TESSELLATION::Generate_Triangles(implicit_surface);surface->Update_Triangle_List();return surface;}

    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {return implicit_surface.box;}

//#####################################################################
};
}
#endif
#endif 
