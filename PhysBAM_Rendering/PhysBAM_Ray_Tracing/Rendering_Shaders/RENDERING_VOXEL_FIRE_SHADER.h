//#####################################################################
// Copyright 2003-2009, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_VOXEL_FIRE_SHADER__
#define __RENDERING_VOXEL_FIRE_SHADER__
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/BLACKBODY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_VOXEL_FIRE_SHADER:public VOLUMETRIC_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    VECTOR<T,3> absorption;
    VECTOR<T,3> scattering;
    T inscattering_amplification_factor;
    bool isotropic_scattering;
    T phase_function_g;
    T emission_amplification;
    T temperature_scale,temperature_offset;
    bool clamp_low_temperature;
    T temperature_lowest;
    T density_scale,density_offset;
    bool clamp_low_density;
    T density_lowest;
    BLACKBODY<T> blackbody;
    MATRIX<T,3> world_xyz_to_display_xyz;
    bool use_lms_scaling;
    LEVELSET_3D<GRID<TV> >* blue_core_levelset;
    RENDER_WORLD<T>& world;
    
    RENDERING_VOXEL_FIRE_SHADER(const VECTOR<T,3>& absorption_input,const VECTOR<T,3>& scattering_input,
        const T inscattering_amplification_factor_input,const T emission_amplification,
        const T temperature_scale_input,const T temperature_offset_input,const bool clamp_low_temperature_input,
        const T temperature_lowest_input,const T white_point_temperature,
        const T density_scale_input,const T density_offset_input,const bool clamp_low_density_input,
        const T density_lowest_input,const bool use_lms_scaling,LEVELSET_3D<GRID<TV> >* blue_core_levelset_input,RENDER_WORLD<T>& world_input) 
        :VOLUMETRIC_SHADER<T>(world_input),absorption(absorption_input),scattering(scattering_input),
         inscattering_amplification_factor(inscattering_amplification_factor_input),isotropic_scattering(false),phase_function_g((T)0.67),
         emission_amplification(emission_amplification),
         temperature_scale(temperature_scale_input),temperature_offset(temperature_offset_input),
         clamp_low_temperature(clamp_low_temperature_input),temperature_lowest(temperature_lowest_input),
         density_scale(density_scale_input),density_offset(density_offset_input),
         clamp_low_density(clamp_low_density_input),density_lowest(density_lowest_input),
         use_lms_scaling(use_lms_scaling),
         blue_core_levelset(blue_core_levelset_input),world(world_input)
    {
        isotropic_scattering=true;
        // init cell min and max on levelset so we can do inside test
        if(blue_core_levelset) blue_core_levelset->Compute_Cell_Minimum_And_Maximum();
        // compute adaption using Von Kries Transform or simple luminance scaling
        VECTOR<T,3> xyz_white_point=blackbody.Calculate_XYZ(white_point_temperature);
        if(use_lms_scaling){
            LOG::cout<<"*****************Using lms"<<std::endl;
            MATRIX<T,3> xyz_to_lms(T(0.3897),T(-0.2298),0,T(.6890),T(1.1834),0,T(-.0787),T(.0464),1);
            MATRIX<T,3> lms_to_xyz=xyz_to_lms.Inverse();
            VECTOR<T,3> lms_scale=xyz_to_lms*xyz_white_point;
            MATRIX<T,3> world_lms_to_display_lms(T(1)/lms_scale.x,0,0,0,T(1)/lms_scale.y,0,0,0,T(1)/lms_scale.z);
            world_xyz_to_display_xyz=lms_to_xyz*world_lms_to_display_lms*xyz_to_lms;}
        else{
            LOG::cout<<"*****************NOT Using lms"<<std::endl;
            world_xyz_to_display_xyz=MATRIX<T,3>(T(1)/xyz_white_point.x,0,0,0,T(1)/xyz_white_point.y,0,0,0,T(1)/xyz_white_point.z);}
        LOG::cout<<"absorption="<<absorption<<", scattering="<<scattering<<", inscattering_amplification_factor="<<inscattering_amplification_factor<<std::endl;
        LOG::cout<<"isotropic_scattering="<<isotropic_scattering<<", phase_function_g="<<phase_function_g<<std::endl;
        LOG::cout<<"emission_amplification="<<emission_amplification<<", temperature_scale="<<temperature_scale<<", temperature_offset="<<temperature_offset<<std::endl;
        LOG::cout<<"clamp_low_temperature="<<clamp_low_temperature<<", temperature_lowest="<<temperature_lowest<<std::endl;
        LOG::cout<<"density_scale="<<density_scale<<", density_offset="<<density_offset<<", clamp_low_density="<<clamp_low_density<<", density_lowest="<<density_lowest<<std::endl;
    }

    virtual ~RENDERING_VOXEL_FIRE_SHADER()
    {}

    T Phase(const VECTOR<T,3>& incoming, const VECTOR<T,3>& outgoing)
    {if(isotropic_scattering)return VOLUMETRIC_SHADER<T>::Isotropic_Phase_Function(incoming,outgoing);
    else return VOLUMETRIC_SHADER<T>::Henyey_Greenstein_Phase_Function(incoming,outgoing,phase_function_g);}
   
    VECTOR<T,3> Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const VECTOR<T,3>& color) PHYSBAM_OVERRIDE
    {const RENDERING_VOXELS<T>* voxel_object=(const RENDERING_VOXELS<T>*)&object;
    T start_t,end_t;
    if(!INTERSECTION::Get_Intersection_Range(ray.ray,voxel_object->box,start_t,end_t))return color;
    T current_t=end_t;T old_t;
    bool last_segment=false;
    VECTOR<T,3> attenuated_color=color;
    int step_count=0;
    while(!last_segment){
        step_count++;
        VECTOR<T,3> current_position=ray.ray.Point(current_t);
        T current_volumetric_step=voxel_object->Volumetric_Integration_Step(RAY<VECTOR<T,3> >(current_position,-ray.ray.direction,true), world.random.Get_Uniform_Number((T)0,(T)1));
        if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Volumetric_Step %f, Current_T %f\n",current_volumetric_step,current_t));
        if(current_t-current_volumetric_step<start_t){last_segment=true;current_volumetric_step=current_t-start_t;}
        VECTOR<T,3> midpoint=ray.ray.Point(current_t-(T)0.5*current_volumetric_step);        
        T density,temperature;
        if(blue_core_levelset && blue_core_levelset->Lazy_Inside(midpoint)){density=0;temperature=289;}else{density=voxel_object->Source_Term(1,midpoint);temperature=voxel_object->Source_Term(2,midpoint);}
        temperature=temperature_scale*(temperature+temperature_offset);
        if(clamp_low_temperature) temperature=max(temperature_lowest,temperature);
        density=density_scale*(density+density_offset);
        if(clamp_low_density) density=max(density_lowest,density);
        // compute emitted radiance
        VECTOR<T,3> emitted_radiance=blackbody.cie.XYZ_To_RGB(world_xyz_to_display_xyz*(blackbody.Calculate_XYZ(temperature)))*emission_amplification;
        // compute extinction
        VECTOR<T,3> attenuation_coefficient=-current_volumetric_step*density*(absorption+scattering);
        VECTOR<T,3> attenuation(exp(attenuation_coefficient.x),exp(attenuation_coefficient.y),exp(attenuation_coefficient.z));
        attenuated_color=attenuated_color*attenuation+density*absorption*emitted_radiance*current_volumetric_step;
        old_t=current_t;current_t-=current_volumetric_step;}
    return attenuated_color;}

    VECTOR<T,3> Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,const VECTOR<T,3>& light_color) PHYSBAM_OVERRIDE
    {return light_color;}

//#####################################################################
};
}
#endif
