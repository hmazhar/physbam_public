//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDER_WORLD 
//#####################################################################
#ifndef __RENDER_WORLD__
#define __RENDER_WORLD__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/CAMERA.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RGB_COLORS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_UNIFORM_COLOR_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
namespace PhysBAM{

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T> class IRRADIANCE_CACHE;
#endif
template<class T> class RENDERING_LIGHT;
template<class T> class RENDERING_OBJECT;
template<class T> class RAY_TRACER_DEBUG_DATA;
template<class T> class VOLUMETRIC_SHADER;
template<class T> class MATERIAL_SHADER;

template<class T>
class RENDER_WORLD:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    RENDERING_OBJECT<T>* ether;
    ARRAY<RENDERING_LIGHT<T>*> lights;
    ARRAY<RENDERING_OBJECT<T>*> standard_objects;
    ARRAY<RENDERING_OBJECT<T>*> volumetric_objects;
    RGB_COLORS<T> rgb_colors;
    CAMERA<T> camera;
    RANDOM_NUMBERS<T> random;
    int ray_depth_limit;
    int photon_depth_limit;
    T ray_contribution_limit;
    bool debug_mode;
    bool use_photon_mapping;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    bool use_irradiance_cache;
#endif
    bool use_adaptive_supersampling;
    int adaptive_supersampling_depth_limit;
    T adaptive_supersampling_tolerance;
    T irradiance_cache_max_distance;
    T max_photon_distance;
    int number_of_photons_for_estimate;
    PHOTON_MAP<T> global_photon_map,caustic_photon_map,volume_photon_map;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    IRRADIANCE_CACHE<T> irradiance_cache;
#endif
    MATERIAL_SHADER<T>* background_shader;
    RENDERING_UNIFORM_COLOR_SHADER<T>* default_background_shader;
    int threads;

    RENDER_WORLD();

    void Clean_Scene()
    {lights.Clean_Memory();standard_objects.Clean_Memory();volumetric_objects.Clean_Memory();standard_objects.Append(ether);}

    void Add_Light(RENDERING_LIGHT<T>* light)
    {light->light_index=lights.m+1;lights.Append(light);}
    
    void Add_Object(RENDERING_OBJECT<T>* object)
    {if(object->material_shader)standard_objects.Append(object);
    if(object->volumetric_shader)volumetric_objects.Append(object);}

    void Set_Ether_Volumetric_Shader(VOLUMETRIC_SHADER<T>* new_shader)
    {VOLUMETRIC_SHADER<T>* old_shader=ether->volumetric_shader;
    ether->volumetric_shader=new_shader;if(!old_shader)volumetric_objects.Append(ether);}

    void Use_Spatial_Partition(const bool use_spatial_partition)
    {}

    void Track_Debugging_Information()
    {debug_mode=true;}

    TV Same_Side_Normal(const RENDERING_RAY<T>& ray,const TV& intersection_point,const RENDERING_OBJECT<T>& object,bool& flipped_normal)
    {TV same_side_normal=object.Normal(intersection_point,ray.ray.aggregate_id);
    if(object.two_sided && TV::Dot_Product(same_side_normal,ray.ray.direction)>0){flipped_normal=true;return -same_side_normal;}else{flipped_normal=false;return same_side_normal;}}

    TV Reflected_Direction(const TV& normal,const TV& incident_ray_direction)
    {return incident_ray_direction-2*TV::Dot_Product(incident_ray_direction,normal)*normal;}

    TV Random_Reflected_Direction_In_Hemisphere(const TV& normal,const T xi_1,const T xi_2)
    {TV pseudo_direction(cos(T(2)*T(pi)*xi_1)*sqrt(xi_2),sin(T(2)*T(pi)*xi_1)*sqrt(xi_2),sqrt(T(1)-xi_2));
    TV temporary_vector=normal.Orthogonal_Vector();
    TV u=TV::Cross_Product(temporary_vector,normal);u.Normalize();
    TV v=TV::Cross_Product(normal,u);
    TV reflected_direction=MATRIX<T,3>(u,v,normal)*pseudo_direction;
    return reflected_direction;}

    const ARRAY<RENDERING_LIGHT<T>*>& Lights()
    {return lights;}

//#####################################################################
    void Render(const RANGE<VECTOR<int,2> >& pixels, PROGRESS_INDICATOR &progress);
    void Render_Pixel(const VECTOR<int,2>& pixel_index,RAY_TRACER_DEBUG_DATA<T>* debug_data=0);
    RENDERING_OBJECT<T>* Closest_Intersection(RENDERING_RAY<T>& ray,const int lowest_priority);
    RENDERING_OBJECT<T>* Point_Inside_Object(const TV& point);
    TV Cast_Ray(RENDERING_RAY<T>& ray,const RENDERING_RAY<T>& parent_ray);
    TV Cast_Ray_For_Photon_Gather(RENDERING_RAY<T>& ray,const RENDERING_RAY<T>& parent_ray);
    TV Incident_Light(RENDERING_RAY<T> ray, const RENDERING_LIGHT<T>& light, const RENDERING_RAY<T>& full_ray,const RENDERING_RAY<T>& parent_ray);
    void Cast_Photon(RENDERING_RAY<T>& ray,const RENDERING_RAY<T>& parent_ray,const TV& power,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,
        const int diffuse_bounces,const int specular_bounces);
    void Use_Photon_Mapping(const int global_photons,const int caustic_photons,const int volume_photons,const T max_photon_distance,const int number_of_photons_for_estimate);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    void Use_Irradiance_Cache(const bool use_irradiance_cache_input,const T irradiance_cache_max_distance_input,const int irradiance_cache_samples_input);
#endif
    void Use_Adaptive_Supersampling(const bool use_adaptive_supersampling_input,const T supersampling_tolerance, const int supersampling_depth_limit);
    void Prepare_For_Forward_Ray_Trace();
    void Shoot_Photons_From_Light_For_Specific_Map(RENDERING_RAY<T>& ray,int light_index,PHOTON_MAP<T>& map,int photon_budget,typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    void Prepare_BSSRDF_Octree();
#endif
 private:
    TV Attenuate_Ray(RENDERING_RAY<T>& ray,const TV& color);
    TV Attenuate_Light_Ray(RENDERING_RAY<T>& ray,const RENDERING_LIGHT<T>& light,const TV& color);
    TV Attenuate_Photon(RENDERING_RAY<T>& ray,const TV& photon_power,bool& should_throw);
//#####################################################################
};
}
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/IRRADIANCE_CACHE.h>
#endif
