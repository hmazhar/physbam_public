//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_VOXELS
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
using namespace PhysBAM;
//#####################################################################
// Function Precompute_Light_Data
//#####################################################################
template<class T> void RENDERING_VOXELS<T>::
Precompute_Light_Data(bool use_fast_precomputation,RENDER_WORLD<T>& world)
{
    LOG::SCOPE scope("precomputation","Precomputing Light Data");
    Prepare_For_Precomputation(world);
    RENDERING_RAY<T> parent_ray;
    ARRAY<VECTOR<T,3> > location_list;
    Get_Node_Locations(location_list);
    ARRAY<int> node_indirection(location_list.m,false);for(int i=1;i<=location_list.m;i++)node_indirection(i)=i;
    ARRAY<T> distance_squared_to_light;if(use_fast_precomputation) distance_squared_to_light.Resize(location_list.m,false,false);
    const ARRAY<RENDERING_LIGHT<T> *>& lights=world.Lights();
    LOG::cout<<"have "<<lights.m<<" lights and "<<location_list.m<<" locations"<<std::endl;
    for(int light_index=1;light_index<=lights.m;light_index++){
        PROGRESS_INDICATOR progress(location_list.m);LOG::cout<<"Light "<<light_index<<": ";
        ARRAY<RAY<VECTOR<T,3> > > sample_array;
        lights(light_index)->Sample_Points(location_list(1),VECTOR<T,3>(1,0,0),sample_array);
        VECTOR<T,3> average_sample_point;for(int sample=1;sample<=sample_array.m;sample++)average_sample_point+=sample_array(sample).Point(sample_array(sample).t_max);
        average_sample_point/=(T)sample_array.m;
        if(use_fast_precomputation){
            for(int i=1;i<=location_list.m;i++)distance_squared_to_light(i)=(average_sample_point-location_list(i)).Magnitude_Squared();
            Sort(node_indirection,Indirect_Comparison(distance_squared_to_light));}
        for(int i=1;i<=location_list.m;i++){
            SEGMENT_3D<T> s(location_list(node_indirection(i)),average_sample_point);
            RAY<VECTOR<T,3> > rr(s);
            RENDERING_RAY<T> ray_to_light(rr,1,this);
            VECTOR<T,3> light_color=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,parent_ray);
            //VECTOR<T,3> light_color=this->volumetric_shader->Attenuate_Light(ray_to_light,*this,*lights(light_index),(*lights(light_index)).Emitted_Light(ray_to_light));
            Set_Precomputed_Light_Data(node_indirection(i),light_index,light_color);
            if(use_fast_precomputation)Set_Precomputed_Light_Valid(node_indirection(i),light_index,true);
            progress.Progress();}}

    Postprocess_Light_Field();

    // Make sure subsequent shading uses precomputation
    if(!use_fast_precomputation)for(int light=1;light<=lights.m;light++)for(int i=1;i<=location_list.m;i++)Set_Precomputed_Light_Valid(i,light,true);
}
//#####################################################################
template class RENDERING_VOXELS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDERING_VOXELS<double>;
#endif
