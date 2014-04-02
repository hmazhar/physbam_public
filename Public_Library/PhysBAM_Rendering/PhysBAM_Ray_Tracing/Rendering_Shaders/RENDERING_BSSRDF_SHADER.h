#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Jiayi Chong, Igor Neverov, Andrew Selle, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// This class is an implementation of Jensen & Buhler, "A Rapid Hierarchical Rendering Technique for Translucent Materials"
// SIGGRAPH 2002.  The variables in the paper map to our variables as follows:
// sigma_tr = effective_transport_coefficient
// sigma_a = absorption coefficient
// sigma_s = scattering coefficient
// sigma_t = extinction coefficient (sigma_a+sigma_t)
// sigma_s' = reduced scattering coefficient
// sigma_t' = reduced extinction coefficient (sigma_s'+sigma_a)
// F_dr = diffuse Fresnel reflectance term
// F_dt = diffuse Fresnel transmittance
// lu = mean free path
// A = boundary_condition_mismatch_correction
// g = phase_function_mean_cosine
// z_r = distance from x0 to dipole light to point of illumination
// z_v = distance from x0 to other dipole light
#ifndef __RENDERING_BSSRDF_SHADER__
#define __RENDERING_BSSRDF_SHADER__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_REFLECTION_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
namespace PhysBAM{

template<class T>
class RENDERING_BSSRDF_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;
  
    const MATERIAL_SHADER<T>* color_shader;
    const MATERIAL_SHADER<T>* ambient_shader;
    T index_of_refraction;
    VECTOR<T,3> effective_transport_coefficient; // sigma_tr
    VECTOR<T,3> diffuse_mean_free_path; // 1/sigma_tr
    T diffuse_fresnel_transmittance;
    T diffuse_fresnel_reflectance;
    T error_criterion;
    T boundary_condition_mismatch;
    T real_light_to_virtual_light_factor;
    int samples_per_octree_cell;
    bool use_irradiance_cache_file;
    std::string irradiance_cache_filename;

    GRID<VECTOR<T,1> > reflectance_grid;
    ARRAY<T,VECTOR<int,1> > reflectance_to_alpha_prime;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T> alpha_interpolation;
    
//#####################################################################
    RENDERING_BSSRDF_SHADER(const T index_of_refraction_input,const T error_criterion_input,const VECTOR<T,3>& diffuse_mean_free_path,const int samples_per_octree_cell,const MATERIAL_SHADER<T>* color_shader_input,const MATERIAL_SHADER<T>* ambient_shader_input,RENDER_WORLD<T>& world_input);
    VECTOR<T,3> Calculate_Surface_Irradiance(const VECTOR<T,3>& point,const VECTOR<T,3>& same_side_normal,RENDER_WORLD<T>& world,const RENDERING_OBJECT<T>& intersection_object);
    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
#endif
