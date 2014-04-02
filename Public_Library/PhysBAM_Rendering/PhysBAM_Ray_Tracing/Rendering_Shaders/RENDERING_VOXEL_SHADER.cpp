//#####################################################################
// Copyright 2003-2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_SHADER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_VOXEL_SHADER<T>::
RENDERING_VOXEL_SHADER(const TV& absorption_input,const TV& scattering_input,const T inscattering_amplification_input,
    const T emission_amplification_input,const T white_point_temperature,const bool use_lms_scaling_input,
    const bool use_blackbody_ramp_input,const bool use_constant_emission_color_input,const TV constant_emission_color_input,
    RENDER_WORLD<T>& world_input)
    :VOLUMETRIC_SHADER<T>(world_input),absorption(absorption_input),absorption_shadow(absorption_input),scattering(scattering_input),
    inscattering_amplification(inscattering_amplification_input),isotropic_scattering(true),phase_function_g((T)0.67),
    emission_amplification(emission_amplification_input),use_lms_scaling(use_lms_scaling_input),empty_implicit_surface(0),
    levelset_multiple_object(0),non_empty_region(0),use_empty_implicit_surface_for_light_attenuation(false),
    use_temperature_remap(false),use_blackbody_ramp(use_blackbody_ramp_input),
    use_constant_emission_color(use_constant_emission_color_input),constant_emission_color(constant_emission_color_input)
{
    supports_photon_mapping=true;
    int ghost_cells=3;
    if(!use_constant_emission_color && emission_amplification_input>0){
        TV xyz_white_point=blackbody.Calculate_XYZ(white_point_temperature);
        if(use_lms_scaling){
            MATRIX<T,3> xyz_to_lms(T(0.3897),T(-0.2298),0,T(.6890),T(1.1834),0,T(-.0787),T(.0464),1);
            MATRIX<T,3> lms_to_xyz=xyz_to_lms.Inverse();
            TV lms_scale=xyz_to_lms*xyz_white_point;
            MATRIX<T,3> world_lms_to_display_lms(T(1)/lms_scale.x,0,0,0,T(1)/lms_scale.y,0,0,0,T(1)/lms_scale.z);
            world_xyz_to_display_xyz=lms_to_xyz*world_lms_to_display_lms*xyz_to_lms;}
        else{world_xyz_to_display_xyz=MATRIX<T,3>(T(1)/xyz_white_point.x,0,0,0,T(1)/xyz_white_point.y,0,0,0,T(1)/xyz_white_point.z);}
        if(use_blackbody_ramp){
            blackbody_ramp_grid.Initialize(500,0,(T)3000);blackbody_ramp.Resize(blackbody_ramp_grid.Domain_Indices(ghost_cells));
            for(int i=1;i<=blackbody_ramp_grid.counts.x;i++) blackbody_ramp(i)=blackbody.cie.XYZ_To_RGB(world_xyz_to_display_xyz*(blackbody.Calculate_XYZ(blackbody_ramp_grid.X(VECTOR<int,1>(i)).x)));
            BOUNDARY_UNIFORM<GRID<VECTOR<T,1> >,TV>().Fill_Ghost_Cells(blackbody_ramp_grid,blackbody_ramp,blackbody_ramp,0,0,ghost_cells);}}

    LOG::cout<<"absorption="<<absorption<<", scattering="<<scattering<<", inscattering_amplification="<<inscattering_amplification<<std::endl;
    LOG::cout<<"emission_amplification="<<emission_amplification<<", white_point_temperature="<<white_point_temperature<<std::endl;
    LOG::cout<<"use_lms_scaling="<<use_lms_scaling<<", use_blackbody_ramp="<<use_blackbody_ramp<<std::endl;
    LOG::cout<<"use_constant_emission_color="<<use_constant_emission_color<<", constant_emission_color="<<constant_emission_color<<std::endl;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RENDERING_VOXEL_SHADER<T>::
~RENDERING_VOXEL_SHADER()
{}
//#####################################################################
// Function Attenuate_Color
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_VOXEL_SHADER<T>::
Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const TV& color)
{
    const RENDERING_VOXELS<T>* voxel_object=(const RENDERING_VOXELS<T>*)&object;
    T start_t,end_t;
    if(!INTERSECTION::Get_Intersection_Range(ray.ray,voxel_object->box,start_t,end_t))return color;
    T current_t=end_t;T old_t;
    bool last_segment=false;
    const ARRAY<RENDERING_LIGHT<T> *>& lights=world.Lights();
    TV attenuated_color=color;
    int step_count=0;
    while(!last_segment){
        step_count++;
        TV current_position=ray.ray.Point(current_t);
        T current_volumetric_step=voxel_object->Volumetric_Integration_Step(RAY<TV>(current_position,-ray.ray.direction,true), world.random.Get_Uniform_Number((T)0,(T)1));
        if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Volumetric_Step %f, Current_T %f\n",current_volumetric_step,current_t));
        if(current_t-current_volumetric_step<start_t){last_segment=true;current_volumetric_step=current_t-start_t;}
        TV sample_point=ray.ray.Point(current_t-(T).5*current_volumetric_step);

        bool empty_multiple_region=false;
        if(levelset_multiple_object) empty_multiple_region=!levelset_multiple_object->Inside_Region_Only(sample_point,non_empty_region);

        if((!empty_implicit_surface || !empty_implicit_surface->Lazy_Inside(sample_point)) && !empty_multiple_region){
            T volumetric_coefficient=voxel_object->Source_Term(1,sample_point);
            // inscattering
            if(inscattering_amplification>0) for(int light_index=1;light_index<=lights.m;light_index++){
                if(world.volume_photon_map.photons.m){
                    TV normal=-ray.ray.direction;normal.Normalize();
                    attenuated_color+=(current_volumetric_step*inscattering_amplification*
                        (world.volume_photon_map.Irradiance_Estimate(sample_point,normal,world.max_photon_distance,world.number_of_photons_for_estimate,ray,PHOTON_MAP<T>::VOLUME_PHOTON_MAP))/(absorption+scattering));}
                TV accumulated_samples(0,0,0);
                if(voxel_object->Use_Precomputed_Light_Data(sample_point,light_index)) {
                    ARRAY<RAY<TV> > sample_array;lights(light_index)->Sample_Points(sample_point,TV(1,0,0),sample_array);
                    accumulated_samples=voxel_object->Precomputed_Light_Data(sample_point,light_index)*Phase(-sample_array(1).direction,-ray.ray.direction);}
                else{
                    ARRAY<RAY<TV> > sample_array;
                    lights(light_index)->Sample_Points(sample_point,TV(1,0,0),sample_array);
                    for(int sample=1;sample<=sample_array.m;sample++){
                        RENDERING_RAY<T> ray_to_light(sample_array(sample),1,&object);
                        TV light_color=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,ray);
                        accumulated_samples+=light_color*Phase(-ray_to_light.ray.direction,-ray.ray.direction);}
                    accumulated_samples/=(T)sample_array.m;}
                attenuated_color+=current_volumetric_step*accumulated_samples*volumetric_coefficient*scattering*inscattering_amplification/T(4*pi);}
            // extinction == absorption and outscattering 
            TV attenuation_coefficient=-current_volumetric_step*volumetric_coefficient*(absorption+scattering);
            TV attenuation(exp(attenuation_coefficient.x),exp(attenuation_coefficient.y),exp(attenuation_coefficient.z));
            attenuated_color*=attenuation;
            // blackbody emission
            if(emission_amplification>0){
                TV emitted_radiance;
                if(use_constant_emission_color) emitted_radiance=constant_emission_color*emission_amplification*volumetric_coefficient*current_volumetric_step;
                else{
                    T temperature=voxel_object->Source_Term(2,sample_point);
                    if(use_temperature_remap)temperature=temperature_remap_interpolation.Clamped_To_Array_Node(temperature_remap_grid,temperature_remap,VECTOR<T,1>(temperature));
                    else if(use_blackbody_ramp){
                        emitted_radiance=blackbody_interpolation.Clamped_To_Array_Node(blackbody_ramp_grid,blackbody_ramp,VECTOR<T,1>(temperature))*emission_amplification*volumetric_coefficient*current_volumetric_step;}
                    else{
                        emitted_radiance=blackbody.cie.XYZ_To_RGB(world_xyz_to_display_xyz*(blackbody.Calculate_XYZ(temperature)))*emission_amplification*volumetric_coefficient*current_volumetric_step;}}
                attenuated_color+=emitted_radiance;}}
        old_t=current_t;current_t-=current_volumetric_step;}
    return attenuated_color;
}
//#####################################################################
// Function Attenuate_Light
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_VOXEL_SHADER<T>::
Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,const TV& light_color)
{
    const RENDERING_VOXELS<T>* voxel_object=(const RENDERING_VOXELS<T>*)&object;
    T start_t,end_t;
    if(!INTERSECTION::Get_Intersection_Range(ray.ray,voxel_object->box,start_t,end_t))return light_color;
    T current_t=start_t;T old_t;
    bool last_segment=false;
    TV transmittance(1,1,1),emitted_radiance(0,0,0);
    int step_count=0;
    while(!last_segment){
        step_count++;
        TV current_position=ray.ray.Point(current_t);
        T current_volumetric_step=voxel_object->Volumetric_Integration_Step(RAY<TV>(current_position,ray.ray.direction,true), world.random.Get_Uniform_Number((T)0,(T) 1));
        if(current_t+current_volumetric_step>end_t){last_segment=true;current_volumetric_step=end_t-current_t;}
        TV sample_point=ray.ray.Point(current_t+(T)0.5*current_volumetric_step);

        bool empty_multiple_region=false;
        if(levelset_multiple_object) empty_multiple_region=!levelset_multiple_object->Inside_Region_Only(sample_point,non_empty_region);

        if((!use_empty_implicit_surface_for_light_attenuation || !empty_implicit_surface || !empty_implicit_surface->Lazy_Inside(sample_point)) && !empty_multiple_region){
            // scattering absorption
            if(ray.current_object->volumetric_shader&&voxel_object->Use_Precomputed_Light_Data(sample_point,light.light_index)){
                T volumetric_coefficient=voxel_object->Source_Term(1,sample_point);
                // absorption and outscattering
                TV attenuation_coefficient=-(T).5*current_volumetric_step*volumetric_coefficient*(absorption_shadow+scattering);
                TV attenuation(exp(attenuation_coefficient.x),exp(attenuation_coefficient.y),exp(attenuation_coefficient.z));
                transmittance*=attenuation;
                // early exit
                TV precomputed_light=voxel_object->Precomputed_Light_Data(sample_point,light.light_index);
                if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Exit early at %f %f %f with precomputed value %f %f %f\n",sample_point.x,sample_point.y,sample_point.z,precomputed_light.x,precomputed_light.y,precomputed_light.z));
                return transmittance*precomputed_light+emitted_radiance;}
            else{
                T volumetric_coefficient=voxel_object->Source_Term(1,sample_point);
                // absorption and outscattering
                TV attenuation_coefficient=-current_volumetric_step*volumetric_coefficient*(absorption_shadow+scattering);
                TV attenuation(exp(attenuation_coefficient.x),exp(attenuation_coefficient.y),exp(attenuation_coefficient.z));
                transmittance*=attenuation;}}
        old_t=current_t;current_t+=current_volumetric_step;}
    //LOG::cout<<"Did "<<step_count<<" steps in volumetric voxels"<<std::endl;
    return transmittance*light_color+emitted_radiance;
}
//#####################################################################
// Function Scatter_Photon_Ray
//#####################################################################
template<class T> bool RENDERING_VOXEL_SHADER<T>::
Scatter_Photon_Ray(const RENDERING_OBJECT<T>& object,RENDERING_RAY<T>& ray,TV& photon_power,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const T fixed_step_size)
{
    assert(!empty_implicit_surface); // not handled for photon mapping yet
    const RENDERING_VOXELS<T>* voxel_object=(const RENDERING_VOXELS<T>*)&object;
    // use russian roulette to determine whether we absorb or scatter
    T source_term=(T).5*(voxel_object->Source_Term(1,ray.ray.endpoint)+voxel_object->Source_Term(1,ray.ray.Point(fixed_step_size)));
    TV extinction_coefficient=source_term*(absorption+scattering);
    T k=(T)one_third*(extinction_coefficient.x+extinction_coefficient.y+extinction_coefficient.z);
    T p=exp(-k*fixed_step_size);
    TV scattering_albedo=scattering/extinction_coefficient;
    T average_scattering_albedo=(T)one_third*(scattering_albedo.x+scattering_albedo.y+scattering_albedo.z);
    T u=world.random.Get_Uniform_Number((T)0,(T)1);
    if(u>p){// scatter
        T p_s=(T)one_third*(scattering/(absorption+scattering)).Sum();
        if(world.random.Get_Uniform_Number((T)0,(T)1)<p_s){
            if(type==PHOTON_MAP<T>::VOLUME_PHOTON_MAP) world.volume_photon_map.Store_Photon(ray.ray.Point((T)-1*log(u)/k),ray.ray.direction,photon_power,ray.recursion_depth);
            photon_power=photon_power*scattering_albedo*(T(1)/T(average_scattering_albedo));
            ray.ray.endpoint=ray.ray.Point(fixed_step_size);
            ray.ray.direction=Inverse_Phase_Direction(ray.ray.direction);
            ray.ray.semi_infinite=true;
            ray.ray.t_max=0;
            return true;}
        else
            return false;}
    else{ray.ray.endpoint=ray.ray.Point(fixed_step_size);ray.ray.semi_infinite=true;return true;}
}
//#####################################################################
template class RENDERING_VOXEL_SHADER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDERING_VOXEL_SHADER<double>;
#endif
