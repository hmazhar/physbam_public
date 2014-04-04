//#####################################################################
// Copyright 2004-2008, Zhaosheng Bao, Eilene Hao, Sergey Koltakov, Frank Losasso, Craig Schroeder, Andrew Selle, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_DIRECTIONAL_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_POINT_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_RECTANGLE_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_SPOTLIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_VOXEL_FIRE_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_LIGHT_SHADER.h>
#include "GENERIC_RENDER_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Light - Prepare lights
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Light(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    int light_index=100;
    std::string type=parameters.Get_Parameter("Type",std::string("Null"));
    bool global_photons=parameters.Get_Parameter("Supports_Global_Photons",true);bool caustic_photons=parameters.Get_Parameter("Supports_Caustic_Photons",true);
    bool volume_photons=parameters.Get_Parameter("Supports_Volume_Photons",true);
    bool photon_source_only=parameters.Get_Parameter("Photon_Source_Only",false);
    bool casts_shadows=parameters.Get_Parameter("Casts_Shadows",true);
    if(type=="Point"){
        TV position=current_transform.Homogeneous_Times(parameters.Get_Parameter("Position",TV()));
        TV color=current_transform.Homogeneous_Times(parameters.Get_Parameter("Color",TV(1,1,1)));
        T power=parameters.Get_Parameter("Power",(T)1);
        lights.push_back(new RENDERING_POINT_LIGHT<T>(position,color,power,true,world,global_photons,caustic_photons,volume_photons));
        LOG::cout<<"Light  Point position="<<position<<" color="<<color<<" power="<<power<<std::endl;}
    else if(type=="Spotlight"){
        TV position=current_transform.Homogeneous_Times(parameters.Get_Parameter("Position",TV()));
        TV color=current_transform.Homogeneous_Times(parameters.Get_Parameter("Color",TV(1,1,1)));
        T power=parameters.Get_Parameter("Power",(T)1);
        TV direction=current_transform.Transposed().Inverse().Homogeneous_Times(parameters.Get_Parameter("Direction",TV(0,0,1)));
        T cone_angle=parameters.Get_Parameter("Cone_Angle",(T)pi/4);
        T penumbra_angle=parameters.Get_Parameter("Penumbra_Angle",(T)pi/4);
        lights.push_back(new RENDERING_SPOTLIGHT<T>(position,color,power,direction,cone_angle,penumbra_angle,world,global_photons,caustic_photons,volume_photons,photon_source_only));}
    else if(type=="Rectangle"){
        TV position=current_transform.Homogeneous_Times(parameters.Get_Parameter("Position",TV()));
        TV color=current_transform.Homogeneous_Times(parameters.Get_Parameter("Color",TV(1,1,1)));
        T power=parameters.Get_Parameter("Power",(T)1);
        TV u_vector=current_transform.Homogeneous_Times(parameters.Get_Parameter("U_Vector",TV(0,0,1)));
        TV v_vector=current_transform.Homogeneous_Times(parameters.Get_Parameter("V_Vector",TV(1,0,0)));
        int u_samples=parameters.Get_Parameter("U_Samples",(int)5);
        int v_samples=parameters.Get_Parameter("V_Samples",(int)5);
        RENDERING_RECTANGLE_LIGHT<T>* light=new RENDERING_RECTANGLE_LIGHT<T>(position,color,power,u_vector,v_vector,u_samples,v_samples,world,global_photons,caustic_photons,volume_photons,true);
        char area_light_identifier[1024];
        lights.push_back(light);
        sprintf(area_light_identifier,"area_light_%d_shader",light_index);
        RENDERING_LIGHT_SHADER<T>* shader=new RENDERING_LIGHT_SHADER<T>(*light,world);
        shaders.Set(area_light_identifier,shader);
        light->material_shader=shader;
        sprintf(area_light_identifier,"area_light_%d",light_index);
        objects.Set(area_light_identifier,light);
        LOG::cout<<"Light  Area  position="<<position<<" color="<<color<<" power="<<power<<std::endl;
        light_index++;}
    else if(type=="Fire_Voxel"){
        std::string voxel_object=parameters.Get_Parameter("Voxel_Object",std::string("<unknown>"));
        std::string fire_shader=parameters.Get_Parameter("Fire_Shader",std::string("<unknown>"));
        RENDERING_OBJECT<T>* rendering_object=0;
        if(!objects.Get(voxel_object,rendering_object)){fprintf(stderr,"Invalid voxel object named '%s' specified\n",voxel_object.c_str());exit(1);}
        VOLUMETRIC_SHADER<T>* volume_shader=0;
        if(!volume_shaders.Get(fire_shader,volume_shader)){fprintf(stderr,"Invalid fire shader named '%s' specified\n",fire_shader.c_str());exit(1);}
        lights.push_back(new RENDERING_VOXEL_FIRE_LIGHT<T>(*(RENDERING_UNIFORM_VOXELS<T>*)rendering_object,
            *(RENDERING_VOXEL_SHADER<T>*)volume_shader,world,global_photons,caustic_photons,volume_photons));}
    else if(type=="Directional") {
        TV direction=current_transform.Homogeneous_Times(parameters.Get_Parameter("Direction",TV()));
        TV color=current_transform.Homogeneous_Times(parameters.Get_Parameter("Color",TV(1,1,1)));
        T power=parameters.Get_Parameter("Power",(T)1);
        lights.push_back(new RENDERING_DIRECTIONAL_LIGHT<T>(direction,color,power,world));
        LOG::cout<<"Light  Directional direction="<<direction<<" color="<<color<<" power="<<power<<std::endl;}
    else PHYSBAM_FATAL_ERROR("Unknown light type:"+type);
    lights.back()->casts_shadows=casts_shadows;
}
//#####################################################################
template void GENERIC_RENDER_EXAMPLE<float,float>::Light(RENDER_WORLD<float>& world,const int frame,PARAMETER_LIST& parameters);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void GENERIC_RENDER_EXAMPLE<double,float>::Light(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,double>::Light(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
#endif
