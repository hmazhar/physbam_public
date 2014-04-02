//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_HOMOGENEOUS_VOLUME_SHADER__
#define __RENDERING_HOMOGENEOUS_VOLUME_SHADER__
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_HOMOGENEOUS_VOLUME_SHADER:public VOLUMETRIC_SHADER<T>
{
public:
    using VOLUMETRIC_SHADER<T>::world;using VOLUMETRIC_SHADER<T>::supports_photon_mapping;

    T absorption,scattering,extinction;
    T volumetric_step;
    
    RENDERING_HOMOGENEOUS_VOLUME_SHADER(const T absorption_input,const T scattering_input,const T volumetric_step_input,RENDER_WORLD<T>& world_input) 
        :VOLUMETRIC_SHADER<T>(world_input),absorption(absorption_input),scattering(scattering_input),volumetric_step(volumetric_step_input)
    {
        extinction=absorption+scattering;
        supports_photon_mapping=false;
    }

    T Phase(const VECTOR<T,3>& incoming, const VECTOR<T,3>& outgoing)
    {return VOLUMETRIC_SHADER<T>::Isotropic_Phase_Function(incoming,outgoing);}

    VECTOR<T,3> Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const VECTOR<T,3>& color) PHYSBAM_OVERRIDE
    {T start_t,end_t;
    if(!object.Get_Intersection_Range(ray.ray,start_t,end_t))return color;
    T current_t=end_t;T old_t;
    bool last_segment=false;
    const ARRAY<RENDERING_LIGHT<T> *>& lights=world.Lights();
    VECTOR<T,3> attenuated_color=color;
    int step_count=0;
    while(!last_segment){
        step_count++;
        T current_volumetric_step=volumetric_step;
        if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Volumetric_Step %f, Current_T %f\n",current_volumetric_step,current_t));
        if(current_t-current_volumetric_step<start_t){last_segment=true;current_volumetric_step=current_t-start_t;}
        VECTOR<T,3> midpoint=ray.ray.Point(current_t-T(.5)*current_volumetric_step);        
        // absorption and outscattering
        T exponential_attenuation=exp(-current_volumetric_step*extinction);
        VECTOR<T,3> attenuation(exponential_attenuation,exponential_attenuation,exponential_attenuation);
        attenuated_color*=attenuation;
        // Compute inscattering (only direct lighting i.e. signle bounce)
        for(int light_index=1;light_index<=lights.m;light_index++){
            VECTOR<T,3> accumulated_radiance_samples(0,0,0);
            ARRAY<RAY<VECTOR<T,3> > > sample_array;
            lights(light_index)->Sample_Points(midpoint,VECTOR<T,3>(1,0,0),sample_array);
            for(int sample=1;sample<=sample_array.m;sample++){
                RENDERING_RAY<T> ray_to_light(sample_array(sample),1,&object);
                VECTOR<T,3> attenuated_light_radiance=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,ray);
                accumulated_radiance_samples+=attenuated_light_radiance*Phase(-ray_to_light.ray.direction,-ray.ray.direction);}
            accumulated_radiance_samples/=(T)sample_array.m;
            attenuated_color+=current_volumetric_step*accumulated_radiance_samples*scattering;}
        // step next
        old_t=current_t;current_t-=current_volumetric_step;}
    return attenuated_color;}

    VECTOR<T,3> Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,const VECTOR<T,3>& light_color) PHYSBAM_OVERRIDE
    {const RENDERING_VOXELS<T>* voxel_object=(const RENDERING_VOXELS<T>*)&object;
    T start_t,end_t;if(!INTERSECTION::Get_Intersection_Range(ray.ray,voxel_object->box,start_t,end_t))return light_color;
    T attenuation_length=end_t-start_t;return std::exp(-attenuation_length*extinction)*light_color;}

    bool Scatter_Photon_Ray(const RENDERING_OBJECT<T>& object,RENDERING_RAY<T>& ray,VECTOR<T,3>& photon_power,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const T fixed_step_size) PHYSBAM_OVERRIDE
    {return false;}

//#####################################################################
};
}
#endif
