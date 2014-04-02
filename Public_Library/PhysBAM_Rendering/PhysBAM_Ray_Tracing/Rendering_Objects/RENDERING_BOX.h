//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX_3D
//#####################################################################
#ifndef __RENDERING_BOX__
#define __RENDERING_BOX__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_BOX:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    BOX<TV> box;

    RENDERING_BOX()
    {}

    RENDERING_BOX(const T xmin,const T xmax,const T ymin,const T ymax,const T zmin,const T zmax)
        :box(xmin,xmax,ymin,ymax,zmin,zmax)
    {}

    RENDERING_BOX(const RANGE<TV>& box)
        :box(box)
    {}

    virtual ~RENDERING_BOX()
    {}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,box,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const PHYSBAM_OVERRIDE
    {if(aggregate == 1){object_tangent=TV(0,0,1);object_bitangent=TV(0,1,0);} // return TV(-1,0,0);
    else if(aggregate == 2){object_tangent=TV(0,1,0);object_bitangent=TV(0,0,1);} // return TV(1,0,0);
    else if(aggregate == 3){object_tangent=TV(1,0,0);object_bitangent=TV(0,0,1);} //  return TV(0,-1,0);
    else if(aggregate == 4){object_tangent=TV(0,0,1);object_bitangent=TV(1,0,0);} //  return TV(0,1,0);
    else if(aggregate == 5){object_tangent=TV(0,1,0);object_bitangent=TV(1,0,0);} //  return TV(0,0,-1);
    else{object_tangent=TV(1,0,0);object_bitangent=TV(0,1,0);}}

    void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const PHYSBAM_OVERRIDE
    {TV object_tangent,object_bitangent;
    Get_Object_Space_Tangent_And_Bitangent(Object_Space_Point(world_space_point),Object_Space_Vector(world_space_normal),aggregate,object_tangent,object_bitangent);
    world_tangent=World_Space_Vector(object_tangent);world_bitangent=World_Space_Vector(object_bitangent);}

    void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const PHYSBAM_OVERRIDE
    {s=T();t=T();}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(box.Normal(aggregate));}

    bool Inside(const TV& location) const  PHYSBAM_OVERRIDE
    {return box.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return box.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return box.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const PHYSBAM_OVERRIDE
    {return World_Space_Point(box.Surface(Object_Space_Point(location)));}

    bool Get_Intersection_Range(const RAY<TV>& ray,T& start_t,T& end_t) const PHYSBAM_OVERRIDE
    {return INTERSECTION::Get_Intersection_Range(Object_Space_Ray(ray),box,start_t,end_t);}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {return TESSELLATION::Generate_Triangles(box);}

//#####################################################################
};
}
#endif
