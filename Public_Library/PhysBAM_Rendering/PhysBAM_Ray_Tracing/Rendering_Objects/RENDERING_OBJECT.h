//#####################################################################
// Copyright 2002-2007, Jiayi Chong, Ronald Fedkiw, Eran Guendelman, Igor Neverov, Andrew Selle, Mike Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_OBJECT
//#####################################################################
#ifndef __RENDERING_OBJECT__
#define __RENDERING_OBJECT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT_ACCELERATION_PRIMITIVE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
#include <iostream>
#include <string>
namespace PhysBAM{

template<class T> class MATERIAL_SHADER;
template<class T> class VOLUMETRIC_SHADER;
template<class T> class RENDERING_BSSRDF_SHADER;
template<class T> class RENDER_WORLD;

template<class T>
class RENDERING_OBJECT:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    T small_number;
    MATRIX<T,4> transform,inverse_transform; // transforms from standard position into actual position
    MATRIX<T,4> solid_texture_transform;
    bool add_to_spatial_partition;
    std::string name;
    T index_of_refraction;
    mutable int priority;
    bool support_transparent_overlapping_objects;
    MATERIAL_SHADER<T>* material_shader;
    VOLUMETRIC_SHADER<T>* volumetric_shader;
    bool two_sided;
    bool flip_normal;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>* bssrdf_tree;
#endif
    RENDERING_BSSRDF_SHADER<T>* bssrdf_shader;

    RENDERING_OBJECT()
        :add_to_spatial_partition(false),index_of_refraction(1),priority(0),support_transparent_overlapping_objects(false),material_shader(0),volumetric_shader(0),two_sided(true),
         flip_normal(false),
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        bssrdf_tree(0),
#endif
        bssrdf_shader(0)
    {
        small_number=Default_Small_Number();
        transform=MATRIX<T,4>::Identity_Matrix();inverse_transform=MATRIX<T,4>::Identity_Matrix();solid_texture_transform=MATRIX<T,4>::Identity_Matrix();
    }

    virtual ~RENDERING_OBJECT()
    {}

    virtual T Index_Of_Refraction(const TV& world_space_location) const
    {return index_of_refraction;}

    virtual void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

    virtual TV Get_Solid_Texture_Coordinates(const TV& object_space_point,const int aggregate) const
    {return solid_texture_transform.Homogeneous_Times(object_space_point);}

    virtual void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

    virtual void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

    void Update_Transform(const MATRIX<T,4>& A)
    {transform=A*transform;inverse_transform=transform.Inverse();}

    void Set_Transform(const MATRIX<T,4>& A)
    {transform=A;inverse_transform=transform.Inverse();}

    RAY<TV> Object_Space_Ray(const RAY<VECTOR<T,3> >& world_space_ray) const
    {RAY<TV> transformed_ray(Object_Space_Point(world_space_ray.endpoint),Object_Space_Vector(world_space_ray.direction));
    transformed_ray.semi_infinite=world_space_ray.semi_infinite;transformed_ray.t_max=world_space_ray.t_max;
    transformed_ray.aggregate_id=world_space_ray.aggregate_id;
    return transformed_ray;}

    TV Object_Space_Point(const TV& world_space_point)  const
    {return inverse_transform.Homogeneous_Times(world_space_point);}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return inverse_transform.Extract_Rotation()*world_space_vector;}

    TV World_Space_Point(const TV& object_space_point) const
    {return transform.Homogeneous_Times(object_space_point);}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return transform.Extract_Rotation()*object_space_vector;}

    RANGE<TV> World_Space_Bounding_Box() const
    {return ORIENTED_BOX<TV>(Object_Space_Bounding_Box(),transform).Axis_Aligned_Bounding_Box();}

    virtual bool Closed_Volume() const // indicates whether Inside/Outside are meaningful for this object
    {return true;}

    static T Default_Small_Number()
    {STATIC_ASSERT((T)false);}

    virtual bool Intersection(RAY<TV>& ray,const int lowest_priority,const RENDERING_OBJECT<T>** intersected_object) const
    {if(priority>=lowest_priority&&Intersection(ray)){*intersected_object=this;return true;}
    else return false;}

    virtual bool Inside(const TV& location,RENDERING_OBJECT<T>** intersected_object) const
    {if(support_transparent_overlapping_objects && Inside(location)){
        *intersected_object=const_cast<RENDERING_OBJECT<T>*>(this);return true;}
    return false;}

//#####################################################################
    virtual bool Intersection(RAY<TV>& ray) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Intersection(RAY<TV>& ray,const int aggregate) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Preprocess_Efficiency_Structures(RENDER_WORLD<T>& world){}
    virtual void Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TV Normal(const TV& location,const int aggregate=0) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Inside(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Lazy_Inside(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Outside(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Lazy_Outside(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Boundary(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TV Surface(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Signed_Distance(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Has_Bounding_Box() const {return false;}
    virtual RANGE<TV> Object_Space_Bounding_Box() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Close_To_Open_Surface(const TV& location,const T threshold_distance) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Get_Intersection_Range(const RAY<TV>& ray,T& start_t,T& end_t) const {return false;}
    virtual void Generate_BSSRDF_Tree(RENDER_WORLD<T>& world){}
    virtual TRIANGULATED_SURFACE<T>* Generate_Triangles()const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();return 0;};
//#####################################################################
};
template<> inline float RENDERING_OBJECT<float>::Default_Small_Number(){return (float)1e-6;}
template<> inline double RENDERING_OBJECT<double>::Default_Small_Number(){return 1e-5;}
}
#endif
