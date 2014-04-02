#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Jiayi Chong, Igor Neverov, Andrew Selle, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/FRESNEL.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BSSRDF_SHADER.h>
using namespace PhysBAM;
//#####################################################################
// Function RENDERING_BSSRDF_SHADER
//#####################################################################
template<class T> RENDERING_BSSRDF_SHADER<T>::
RENDERING_BSSRDF_SHADER(const T index_of_refraction_input,const T error_criterion_input,const TV& diffuse_mean_free_path,const int samples_per_octree_cell_input,const MATERIAL_SHADER<T>* color_shader_input,const MATERIAL_SHADER<T>* ambient_shader_input,RENDER_WORLD<T>& world_input) 
    :MATERIAL_SHADER<T>(world_input),color_shader(color_shader_input),ambient_shader(ambient_shader_input),index_of_refraction(index_of_refraction_input),error_criterion(error_criterion_input),samples_per_octree_cell(samples_per_octree_cell_input),use_irradiance_cache_file(false)
{
    diffuse_fresnel_reflectance=(T)-1.440/sqr(index_of_refraction)+(T).710/index_of_refraction+(T).668+(T).0636*index_of_refraction;
    boundary_condition_mismatch=(1+diffuse_fresnel_reflectance)/(1-diffuse_fresnel_reflectance);
    diffuse_fresnel_transmittance=1-diffuse_fresnel_reflectance;
    real_light_to_virtual_light_factor=(1+(T)4/3*boundary_condition_mismatch);
    effective_transport_coefficient=Inverse(diffuse_mean_free_path);
    // build lookup table for reduced albedo (alpha' = sigma_s' / sigma_t')
    INTERPOLATION_CURVE<T,T> curve;GRID<VECTOR<T,1> > alpha_prime_grid(100,0,1);
    for(int i=1;i<=alpha_prime_grid.counts.x;i++){
        T alpha_prime=alpha_prime_grid.Axis_X(i,1);
        T Rd=alpha_prime*(1+exp(-4/3*boundary_condition_mismatch*sqrt(3*(1-alpha_prime))))*exp(-sqrt(3*(1-alpha_prime)));
        curve.Add_Control_Point(Rd,alpha_prime);}
    // create the map from reflectance to alpha_prime
    reflectance_grid.Initialize(1000,0,1);reflectance_to_alpha_prime.Resize(reflectance_grid.Domain_Indices(3),false);
    for(int i=1;i<=reflectance_grid.counts.x;i++) reflectance_to_alpha_prime(i)=curve.Value(reflectance_grid.Axis_X(i,1));
    BOUNDARY_UNIFORM<GRID<VECTOR<T,1> >,T>().Fill_Ghost_Cells(reflectance_grid,reflectance_to_alpha_prime,reflectance_to_alpha_prime,0,0);
}
//#####################################################################
// Function Calculate_Surface_Irradiance
//#####################################################################
template<class T> TV RENDERING_BSSRDF_SHADER<T>::
Calculate_Surface_Irradiance(const TV& point,const TV& same_side_normal,RENDER_WORLD<T>& world,const RENDERING_OBJECT<T>& intersection_object) 
{
    const ARRAY<RENDERING_LIGHT<T>*>& lights=world.Lights();
    TV same_side_position=point+same_side_normal*intersection_object.small_number*2;
    TV accumulated_color(0,0,0);
    RENDERING_RAY<T> dummy_root;dummy_root.recursion_depth=0;dummy_root.current_object=world.ether;
    for(int light_index=1;light_index<=lights.m;light_index++){
        TV accumulated_samples(0,0,0);ARRAY<RAY<TV> > sample_array;lights(light_index)->Sample_Points(same_side_position,same_side_normal,sample_array);
        for(int sample=1;sample<=sample_array.m;sample++){
            RENDERING_RAY<T> ray_to_light(sample_array(sample),1,&intersection_object);
            TV light_color=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,dummy_root);
            T L_N=TV::Dot_Product(ray_to_light.ray.direction,same_side_normal);
            if(L_N<0) continue;
            accumulated_samples+=L_N*light_color;}
        accumulated_color+=accumulated_samples/T(sample_array.m);}
    return accumulated_color;
}
//#####################################################################
// Function Shade_Surface_Using_Direct_Illumination
//#####################################################################
template<class T> TV RENDERING_BSSRDF_SHADER<T>::
Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const
{
    assert(intersection_object.bssrdf_shader&&intersection_object.bssrdf_tree);
    ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> > samples_to_use;samples_to_use.Preallocate(4096);
    intersection_object.bssrdf_tree->Find_Irradiance_Candidates(samples_to_use,intersection_point);
    TV diffuse_color(color_shader->Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal));
    VECTOR<double,3> alpha_prime(alpha_interpolation.Clamped_To_Array(reflectance_grid,reflectance_to_alpha_prime,VECTOR<T,1>(diffuse_color.x)),
        alpha_interpolation.Clamped_To_Array(reflectance_grid,reflectance_to_alpha_prime,VECTOR<T,1>(diffuse_color.y)),
        alpha_interpolation.Clamped_To_Array(reflectance_grid,reflectance_to_alpha_prime,VECTOR<T,1>(diffuse_color.z)));
    VECTOR<double,3> real_light_depth=sqrt((double)3*(VECTOR<double,3>(1,1,1)-alpha_prime))/VECTOR<double,3>(effective_transport_coefficient);
    VECTOR<double,3> virtual_light_depth=real_light_depth*(1+4/3*boundary_condition_mismatch);
    VECTOR<double,3> real_light_depth_squared=sqr(real_light_depth),virtual_light_depth_squared=sqr(virtual_light_depth);


    VECTOR<double,3> radiant_exitance_sum(0,0,0);
    VECTOR<double,3> ambient_irradiance(0,0,0);
    if(ambient_shader) ambient_irradiance=VECTOR<double,3>(ambient_shader->Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal));
    if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("diffuse_color %f %f %f alpha_prime %f %f %f\n",diffuse_color.x,diffuse_color.y,diffuse_color.z,alpha_prime.x,alpha_prime.y,alpha_prime.z));
    for(int i=1;i<=samples_to_use.m;++i){
        if(samples_to_use(i).transmitted_irradiance==TV(0,0,0))continue;
        T sample_distance_squared=(intersection_point-samples_to_use(i).position).Magnitude_Squared();
        VECTOR<double,3> vector_sample_distance_squared(sample_distance_squared,sample_distance_squared,sample_distance_squared);
        VECTOR<double,3> distance_to_real=sqrt(vector_sample_distance_squared+real_light_depth_squared);
        VECTOR<double,3> distance_to_virtual=sqrt(vector_sample_distance_squared+virtual_light_depth_squared);
        VECTOR<double,3> real_attenuation_length=distance_to_real*VECTOR<double,3>(effective_transport_coefficient),virtual_attenuation_length=distance_to_virtual*VECTOR<double,3>(effective_transport_coefficient);
        VECTOR<double,3> dE_dPhi=real_light_depth*(real_attenuation_length+VECTOR<double,3>(1,1,1))*exp(-real_attenuation_length)*Inverse(cube(distance_to_real))
                   +virtual_light_depth*(virtual_attenuation_length+VECTOR<double,3>(1,1,1))*exp(-virtual_attenuation_length)*Inverse(cube(distance_to_virtual));
        radiant_exitance_sum+=(double)(diffuse_fresnel_transmittance)*dE_dPhi*(VECTOR<double,3>(samples_to_use(i).transmitted_irradiance_product)+(double)samples_to_use(i).area*ambient_irradiance);}
    VECTOR<double,3> bssrdf_reflection=(double)diffuse_fresnel_transmittance*radiant_exitance_sum/(double)pi*one_over_four_pi;
    return TV(bssrdf_reflection);
}
//#####################################################################
template class RENDERING_BSSRDF_SHADER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDERING_BSSRDF_SHADER<double>;
#endif
#endif
