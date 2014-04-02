//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Sergey Koltakov, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_PLANE
//#####################################################################
#ifndef __RENDERING_PLANE__
#define __RENDERING_PLANE__

#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_PLANE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    PLANE<T> plane;
    TV texture_vector1,texture_vector2;

    RENDERING_PLANE()
    {}

    RENDERING_PLANE(const TV& normal_input,const TV& x1_input)
        :plane(normal_input,x1_input)
    {}
    
    RENDERING_PLANE(const TV& x1_input,const TV& x2_input,const TV& x3_input)
        :plane(x1_input,x2_input,x3_input)
    {}

    virtual ~RENDERING_PLANE()
    {}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,plane,small_number)){ray.semi_infinite=false;
        ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const PHYSBAM_OVERRIDE
    {object_tangent=texture_vector1;object_bitangent=texture_vector2;}

    void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const PHYSBAM_OVERRIDE
    {world_tangent=World_Space_Vector(texture_vector1);world_bitangent=World_Space_Vector(texture_vector2);}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(plane.Normal());}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return plane.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return plane.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return plane.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const PHYSBAM_OVERRIDE
    {return World_Space_Point(plane.Surface(Object_Space_Point(location)));}

    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {return plane.Signed_Distance(Object_Space_Point(location));}

    void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const PHYSBAM_OVERRIDE
    {TV p=object_space_point-plane.x1;
    s=TV::Dot_Product(p,texture_vector1);t=TV::Dot_Product(p,texture_vector2);}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {TRIANGULATED_SURFACE<T>* surface;
    TV u_vector=texture_vector1,v_vector=texture_vector2;
    GEOMETRY_PARTICLES<TV>* particles=new GEOMETRY_PARTICLES<TV>();
    int vertex_1=particles->array_collection->Add_Element(),vertex_2=particles->array_collection->Add_Element(),vertex_3=particles->array_collection->Add_Element(),vertex_4=particles->array_collection->Add_Element();
    particles->X(vertex_1)=(plane.x1);particles->X(vertex_2)=(plane.x1+u_vector);
    particles->X(vertex_3)=(plane.x1+u_vector+v_vector);particles->X(vertex_4)=(plane.x1+v_vector);
    ARRAY<VECTOR<int,3> > triangles(2);triangles(1).Set(1,2,4);triangles(2).Set(4,2,3);
    TRIANGLE_MESH* mesh=new TRIANGLE_MESH(particles->array_collection->Size(),triangles);
    surface=new TRIANGULATED_SURFACE<T>(*mesh,*particles);
    surface->Update_Triangle_List();surface->Update_Vertex_Normals();return surface;};
    
//#####################################################################
};   
}
#endif

