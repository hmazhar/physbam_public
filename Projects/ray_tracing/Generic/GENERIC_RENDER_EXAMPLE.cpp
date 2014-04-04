//#####################################################################
// Copyright 2004-2007, Zhaosheng Bao, Eilene Hao, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Andrew Selle, Tamar Shinar, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Interpolation/BSPLINE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Read_Write/Rendering_Shaders/READ_WRITE_SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include "GENERIC_RENDER_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Initialize_Scene
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Initialize_Scene(RENDER_WORLD<T>& world,const int frame)
{
    // clear out transformation stack and parser
    current_transform=MATRIX<T,4>::Identity_Matrix();transform_stack.Remove_All();
    parser->Reset_Parser();
    // go through each block
    while(1){
        std::string identifier;PARAMETER_LIST parameters;
        if(!parser->Get_Statement(identifier,parameters))break;
        if(identifier=="Camera")Camera(world,frame,parameters);
        else if(identifier=="Options")Options(world,frame,parameters);
        else if(identifier=="Light")Light(world,frame,parameters);
        else if(identifier=="Material")Material(world,frame,parameters);
        else if(identifier=="Volume_Material")Volume_Material(world,frame,parameters);
        else if(identifier=="Transform")Transform(world,frame,parameters);
        else if(identifier=="Object")Object(world,frame,parameters);
        else if(identifier=="List_Object")List_Object(world,frame,parameters);
        else LOG::cerr<<"Unknown block '"<<identifier<<"'parsed."<<std::endl;}
    // Compute accweleration structures
    List_Object_Compute_Acceleration_Structures();
    // Finally add all
    for(typename std::list<RENDERING_LIGHT<T>*>::iterator i=lights.begin();i!=lights.end();++i){
        LOG::cout<<"Adding Light"<<std::endl;;world.Add_Light((*i));}
    for(typename HASHTABLE<std::string,RENDERING_OBJECT<T>*>::ITERATOR i(objects);i.Valid();i.Next()){
        if(use_spatial_partition && i.Data()->add_to_spatial_partition){
            LOG::cout<<"Adding Object to Spatial Partition "<<i.Key()<<std::endl;
            accelerator.Add_Object(i.Data());}
        else{
            LOG::cout<<"Adding object to normal list "<<i.Key()<<std::endl;
            world.Add_Object(i.Data());}}
    LOG::cout<<"Adding Acceleration Structure..."<<std::endl;
    world.Add_Object(&accelerator);
}
//#####################################################################
// Camera - Prepare render camera and film
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Camera(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    LOG::cout<<"parsing camera"<<std::endl;
    T field_of_view=T(pi)/T(180)*parameters.Get_Parameter("Field_Of_View",T(60));
    T focal_distance=parameters.Get_Parameter("Focal_Distance",T(1.0));
    T aspect_ratio=parameters.Get_Parameter("Aspect_Ratio",T(4.0/3.0));
    gamma_correction=parameters.Get_Parameter("Gamma_Correction",(T)2.2);
    int pixels_x_direction=parameters.Get_Parameter("Width",640);
    int pixels_y_direction=parameters.Get_Parameter("Height",480);
    output_filename=parameters.Get_Parameter("Output_Filename",std::string( "image.%d.rgb"));
    alpha_filename=parameters.Get_Parameter("Alpha_Filename",std::string("alpha.%d.pbm"));
    keep_old_files=parameters.Get_Parameter("Keep_Old_Files",(bool)false);
    // set up camera position, aim and focus
    TV location,look_at,pseudo_up;
    if(parameters.Is_Defined("Motion_Filename")) Get_Camera_Frame(frame,parameters.template Get_Parameter<std::string>("Motion_Filename"),location,look_at,pseudo_up);
    else{
        location=parameters.Get_Parameter("Location",TV(0,0,0));
        look_at=parameters.Get_Parameter("Look_At",TV(0,0,0));
        pseudo_up=parameters.Get_Parameter("Pseudo_Up",TV(0,1,0));}
    PHYSBAM_DEBUG_PRINT("Camera Position/Aim",location,look_at,pseudo_up);
    world.camera.Position_And_Aim_Camera(location,look_at,pseudo_up);
    world.camera.Focus_Camera(focal_distance,aspect_ratio,field_of_view);
    world.camera.film.Set_Resolution(pixels_x_direction,pixels_y_direction);
    PHYSBAM_DEBUG_PRINT("Camera ",world.camera.position,world.camera.focal_point,world.camera.look_vector,world.camera.vertical_vector,world.camera.horizontal_vector);
    // Pixel filtering
    typename FILM<T>::PIXEL_FILTER filter;std::string filter_string=parameters.Get_Parameter("Pixel_Filter",std::string("Gaussian"));
    VECTOR<T,2> filter_width=parameters.Get_Parameter("Pixel_Filter_Width",VECTOR<T,2>(2,2));
    if(filter_string=="Gaussian") filter=Film_Gaussian_Filter<T>;
    else if(filter_string=="Box") filter=Film_Box_Filter<T>;
    else if(filter_string=="Lanczos") filter=Film_Lanczos_Filter<T>;
    else{LOG::cerr<<"Invalid pixel filter"<<std::endl;PHYSBAM_FATAL_ERROR();}
    world.camera.film.Set_Filter(filter_width,filter);
    // Jittering to reduce mach banding from quantization, should probably be lower if using HDR formats
    world.camera.film.dither_amplitude=parameters.Get_Parameter("Dither_Amplitude",(T).5);
    // get clipping region
    clipping_region.min_corner.x=parameters.Get_Parameter("imin",-10000);
    clipping_region.max_corner.x=parameters.Get_Parameter("imax",10000);
    clipping_region.min_corner.y=parameters.Get_Parameter("jmin",-10000);
    clipping_region.max_corner.y=parameters.Get_Parameter("jmax",10000);
}
//#####################################################################
// Camera - Prepare render camera and film
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Get_Camera_Frame(const int camera_frame,const std::string& motion_filename,TV& current_location,TV& current_look_at,TV& current_pseudo_up)
{
    GENERIC_PARSER<T> parser(motion_filename,frame);
    ARRAY<T> frames;

    ARRAY<TV> look_at_array;
    ARRAY<TV> location_array;
    ARRAY<TV> pseudo_up_array;
    while(1){
        std::string identifier;PARAMETER_LIST parameters;
        if(!parser.Get_Statement(identifier,parameters)) break;
        if(identifier!="Camera_Frame")LOG::cerr<<"Get_Camera_Frame(): Unknown block '"<<identifier<<"'parsed."<<std::endl;
        if(!parameters.Is_Defined("Frame")) LOG::cerr<<"Get_Camera_Frame(): specify frame numbers in "<<motion_filename<<std::endl;
        if(!parameters.Is_Defined("Location")) LOG::cerr<<"Get_Camera_Frame(): specify Location in "<<motion_filename<<std::endl;
        if(!parameters.Is_Defined("Look_At")) LOG::cerr<<"Get_Camera_Frame(): specify Look_At in "<<motion_filename<<std::endl;
        if(!parameters.Is_Defined("Pseudo_Up")) LOG::cerr<<"Get_Camera_Frame(): specify Pseudo_Up in "<<motion_filename<<std::endl;
        int keyframe=parameters.template Get_Parameter<int>("Frame");
        LOG::cout<<"Reading camera keyframe: "<<keyframe<<std::endl;
        TV frame_location=parameters.template Get_Parameter<TV>("Location");
        TV frame_look_at=parameters.template Get_Parameter<TV>("Look_At");
        TV frame_pseudo_up=parameters.template Get_Parameter<TV>("Pseudo_Up");
        look_at_array.Append(frame_look_at);
        location_array.Append(frame_location);
        pseudo_up_array.Append(frame_pseudo_up);
        frames.Append((T)keyframe);}
    BSPLINE<T,TV> look_at_spline(frames,look_at_array,3);look_at_spline.Clamp_End_Points();
    BSPLINE<T,TV> location_spline(frames,location_array,3);location_spline.Clamp_End_Points();
    BSPLINE<T,TV> pseudo_up_spline(frames,pseudo_up_array,3);pseudo_up_spline.Clamp_End_Points();
    current_location=location_spline.Clamped_Evaluate((T)camera_frame);
    LOG::cout<<"location"<<current_location<<std::endl;
    current_look_at=look_at_spline.Clamped_Evaluate((T)camera_frame);
    LOG::cout<<"look_at"<<current_look_at<<std::endl;
    current_pseudo_up=pseudo_up_spline.Clamped_Evaluate((T)camera_frame);
    LOG::cout<<"pseudo_up"<<current_pseudo_up<<std::endl;
}
//#####################################################################
// Camera - Prepare render camera and film
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Options(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    world.threads=parameters.Get_Parameter("Threads",(int)1);
    if(world.threads>1) LOG::cout<<"Running multithreaded with "<<world.threads<<" threads"<<std::endl;
    bool use_four_subpixels=parameters.Get_Parameter("Use_Four_Subpixels",(bool)false);
    bool use_jitter=parameters.Get_Parameter("Use_Jitter",(bool)false);
    bool best_candidate=parameters.Get_Parameter("Best_Candidate",(bool)false);
    int samples_per_pixel=parameters.Get_Parameter("Samples_Per_Pixel",(int)1);
    bool use_adaptive_supersampling=parameters.Get_Parameter("Use_Adaptive_Supersampling",(bool)false);
    int adaptive_supersampling_depth_limit=parameters.Get_Parameter("Adaptive_Supersampling_Depth_Limit",(int)20);
    float adaptive_supersampling_tolerance=parameters.Get_Parameter("Adaptive_Supersampling_Tolerance",(float)1e-4);
    world.Use_Adaptive_Supersampling(use_adaptive_supersampling,adaptive_supersampling_tolerance,adaptive_supersampling_depth_limit);
    world.ray_depth_limit=parameters.Get_Parameter("Ray_Depth_Limit",(int)13);
    world.photon_depth_limit=parameters.Get_Parameter("Photon_Depth_Limit",(int)13);
    bool high_quality=parameters.Get_Parameter("High_Quality",(bool)false);
    FILM<T>& film=world.camera.film;
    if(high_quality){samples_per_pixel=4;use_four_subpixels=true;use_jitter=true;}
    film.samples_per_pixel=samples_per_pixel;film.use_four_subpixels=use_four_subpixels;
    if(best_candidate) film.sampler=film.BEST_CANDIDATE;
    else if(use_jitter) film.sampler=film.JITTERED;
    else film.sampler=film.UNIFORM;
    LOG::cout<<"Quality parameters: "<<(use_four_subpixels?"4 subpixels, ":"")<<(use_jitter?"jittered, ":"")<<samples_per_pixel<<" samples per pixel"<<std::endl;
    use_spatial_partition=parameters.Get_Parameter("Spatial_Partition",(bool)true);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(parameters.Get_Parameter("Use_Photon_Map",(bool)false)){
        bool use_irradiance_cache=parameters.Get_Parameter("Use_Irradiance_Cache",(bool)true);
        int global_photons=parameters.Get_Parameter("Global_Photons",(int)0);
        int caustic_photons=parameters.Get_Parameter("Caustic_Photons",(int)0);
        int volume_photons=parameters.Get_Parameter("Volume_Photons",(int)0);
        T max_photon_distance=parameters.Get_Parameter("Max_Photon_Distance",(T)0.15);
        T max_irradiance_cache_distance=parameters.Get_Parameter("Max_Irradiance_Cache_Distance",(T)0.15);
        int irradiance_cache_samples=parameters.Get_Parameter("Irradiance_Cache_Samples",(int)20);
        int number_of_photons_used=parameters.Get_Parameter("Number_Of_Photons_In_Estimate",50);
        world.Use_Photon_Mapping(global_photons,caustic_photons,volume_photons,max_photon_distance,number_of_photons_used);
        world.Use_Irradiance_Cache(use_irradiance_cache,max_irradiance_cache_distance,irradiance_cache_samples);}
#endif
    std::string background_shader=parameters.Get_Parameter("Background_Shader",std::string("<unknown>"));
    if(background_shader!="<unknown>"){
        if(!shaders.Get(background_shader,world.background_shader)){LOG::cerr<<"Invalid shader '"<<background_shader<<"' specified for background"<<std::endl;exit(1);}}
}
//#####################################################################
// Transform - Handle a transformation entry
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Transform(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    std::string type=parameters.Get_Parameter("Type",std::string("Null"));
    if(type!="Null"){
        if(type=="Push")transform_stack.Push(current_transform);
        else if(type=="Pop"){
            if(transform_stack.Empty()){LOG::cerr<<"Transform popped empty stack"<<std::endl;exit(1);}
            current_transform=transform_stack.Pop();}
        else if(type=="Translate"){
            TV vector=parameters.Get_Parameter("Vector",TV(0,0,0));
            current_transform*=MATRIX<T,4>::Translation_Matrix(vector);
            LOG::cout<<"Translating..."<<vector<<std::endl;}
        else if(type=="Rotate"){
            TV axis=parameters.Get_Parameter("Axis",TV(1,0,0));
            T radians;
            if(parameters.Is_Defined("Radians")) radians=parameters.Get_Parameter("Radians",T(0));
            else if(parameters.Is_Defined("Degrees")) radians=(T)pi*parameters.Get_Parameter("Degrees",T(0))/(T)180;
            else{std::cout<<"Failed to get degrees of radians from transform"<<std::endl;exit(1);}
            current_transform*=MATRIX<T,4>::Rotation_Matrix(axis,radians);}
        else if(type=="Scale"){
            TV scale=parameters.Get_Parameter("Vector",TV(1,1,1));
            current_transform*=MATRIX<T,4>::Scale_Matrix(scale);}
        else if(type=="Matrix"){
            TV u=parameters.Get_Parameter("u",TV(1,0,0)),v=parameters.Get_Parameter("v",TV(0,1,0)),
                w=parameters.Get_Parameter("w",TV(0,0,1)),x=parameters.Get_Parameter("x",TV(0,0,0));
            current_transform*=MATRIX<T,4>(u.x,u.y,u.z,0,v.x,v.y,v.z,0,w.x,w.y,w.z,0,x.x,x.y,x.z,1);}}
}
//#####################################################################
// Get_New_Render
//#####################################################################
template<class T,class RW> RENDERING_TRIANGULATED_SURFACE<T>* GENERIC_RENDER_EXAMPLE<T,RW>::
Get_New_Render_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface,const bool smooth_normals,const bool preserve_creases){
    if(smooth_normals) triangulated_surface.Use_Vertex_Normals();
    if(preserve_creases) triangulated_surface.avoid_normal_interpolation_across_sharp_edges=true;
    return new RENDERING_TRIANGULATED_SURFACE<T>(triangulated_surface);
}
//#####################################################################
// Function Add_Solid_Texture
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Add_Solid_Texture(RENDERING_OBJECT<T>* object,PARAMETER_LIST& parameters)
{
    MATRIX<T,4> volume_texture_coordinate=MATRIX<T,4>::Identity_Matrix();
    TV solid_texture_angles=parameters.Get_Parameter("Solid_Texture_Angles",TV());
    TV solid_texture_scale=parameters.Get_Parameter("Solid_Texture_Scale",TV(1,1,1));
    TV solid_texture_translation=parameters.Get_Parameter("Solid_Texture_Translation",TV());
    volume_texture_coordinate*=MATRIX<T,4>::Translation_Matrix(solid_texture_translation);
    volume_texture_coordinate*=MATRIX<T,4>::Rotation_Matrix_X_Axis(solid_texture_angles.x);
    volume_texture_coordinate*=MATRIX<T,4>::Rotation_Matrix_Y_Axis(solid_texture_angles.y);
    volume_texture_coordinate*=MATRIX<T,4>::Rotation_Matrix_Z_Axis(solid_texture_angles.z);
    volume_texture_coordinate*=MATRIX<T,4>::Scale_Matrix(solid_texture_scale);
    LOG::cout<<"Solid Texture angles="<<solid_texture_angles<<" scale="<<solid_texture_scale<<" translation="<<solid_texture_translation<<std::endl;
    object->solid_texture_transform=volume_texture_coordinate;
}
//#####################################################################
#define INSTANTIATE_HELPER(T,RW) \
    template void GENERIC_RENDER_EXAMPLE<T,RW>::Initialize_Scene(RENDER_WORLD<T>& world,const int frame); \
    template void GENERIC_RENDER_EXAMPLE<T,RW>::Camera(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters); \
    template void GENERIC_RENDER_EXAMPLE<T,RW>::Options(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters); \
    template void GENERIC_RENDER_EXAMPLE<T,RW>::Transform(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters); \
    template RENDERING_TRIANGULATED_SURFACE<T>* GENERIC_RENDER_EXAMPLE<T,RW>::Get_New_Render_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface,const bool smooth_normals,const bool preserve_creases); \
    template void GENERIC_RENDER_EXAMPLE<T,RW>::Add_Solid_Texture(RENDERING_OBJECT<T>* object,PARAMETER_LIST& parameters);

INSTANTIATE_HELPER(float,float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATE_HELPER(double,float)
INSTANTIATE_HELPER(double,double)
#endif
