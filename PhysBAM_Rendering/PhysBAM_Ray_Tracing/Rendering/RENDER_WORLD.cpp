//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDER_WORLD
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RAY_TRACER_DEBUG_DATA.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_ETHER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
using namespace PhysBAM;

//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDER_WORLD<T>::
RENDER_WORLD()
    :ray_depth_limit(20),photon_depth_limit(10),ray_contribution_limit((T).001),debug_mode(false),
    use_photon_mapping(false),
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    use_irradiance_cache(false),
#endif
    global_photon_map(0),caustic_photon_map(0),volume_photon_map(0),threads(1)
{
    default_background_shader=new RENDERING_UNIFORM_COLOR_SHADER<T>(RGB_COLORS<T>::Black(),*this);background_shader=default_background_shader;
    ether=new RENDERING_ETHER<T>();
    Clean_Scene();
    random.Set_Seed(23123);
}
//#####################################################################
// Render
//#####################################################################
#ifdef USE_PTHREADS
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
template<class T>
class RENDER_TASK:public THREAD_QUEUE::TASK
{
    RENDER_WORLD<T>& world;
    PROGRESS_INDICATOR& progress;
    pthread_mutex_t* mutex;
    RANGE<VECTOR<int,2> > box;

public:
    RENDER_TASK(RENDER_WORLD<T>& world,PROGRESS_INDICATOR& progress,pthread_mutex_t* mutex,const RANGE<VECTOR<int,2> >& box)
        :world(world),progress(progress),mutex(mutex),box(box)
    {}

    void Run()
    {RENDERING_RAY<T> dummy_root;
    ARRAY<typename FILM<T>::SAMPLE> samples;
    world.camera.film.Generate_Samples(box,world.camera,samples);
    for(int i=1;i<=samples.m;i++){
        typename FILM<T>::SAMPLE& sample=samples(i);
        RENDERING_RAY<T> ray(RAY<VECTOR<T,3> >(world.camera.position,sample.world_position-world.camera.position),1,world.Point_Inside_Object(world.camera.position));
        dummy_root.recursion_depth=0;dummy_root.current_object=world.ether;
        sample.radiance=world.Cast_Ray(ray,dummy_root);
        
        world.camera.film.Add_Sample(sample);
        if(i==samples.m){
            progress.Progress();pthread_mutex_lock(mutex);pthread_mutex_unlock(mutex);}
    }}

    void Commit()
    {progress.Progress();}
};
#endif
template<class T> void RENDER_WORLD<T>::
Render(const RANGE<VECTOR<int,2> >& pixels, PROGRESS_INDICATOR &progress)
{
    RENDERING_RAY<T> dummy_root;
    VECTOR<int,2> box_size=camera.film.Samples_Extent();
    VECTOR<int,2> lengths=pixels.Edge_Lengths();
    ARRAY<typename FILM<T>::SAMPLE> samples;
    progress.Initialize((lengths.x/box_size.x+1)*(lengths.y/box_size.y+1)); // +1 to account for rounding; occasionally will be incorrect, but will get close enough to 100% to judge
                                                                            // progress
    if(threads>1){
#ifdef USE_PTHREADS
        pthread_mutex_t mutex;
        pthread_mutex_init(&mutex,0);
        THREAD_QUEUE task_queue(4);
        for(int box_x=0;box_x*box_size.x<lengths.x;box_x++) for(int box_y=0;box_y*box_size.y<lengths.y;box_y++){
                VECTOR<int,2> lower_corner=VECTOR<int,2>(box_x*box_size.x,box_y*box_size.y)+pixels.min_corner;
                task_queue.Queue(new RENDER_TASK<T>(*this,progress,&mutex,RANGE<VECTOR<int,2> >(lower_corner,lower_corner+box_size)));}
        task_queue.Wait();
        pthread_mutex_destroy(&mutex);
#else
        PHYSBAM_FATAL_ERROR("Threads non 1, but USE_PTHREADS disabled");
#endif
    }
    else for(int box_x=0;box_x*box_size.x<lengths.x;box_x++) for(int box_y=0;box_y*box_size.y<lengths.y;box_y++){
        samples.Remove_All();
        VECTOR<int,2> lower_corner=VECTOR<int,2>(box_x*box_size.x,box_y*box_size.y)+pixels.min_corner;
        camera.film.Generate_Samples(RANGE<VECTOR<int,2> >(lower_corner,lower_corner+box_size),camera,samples);
        progress.Progress();
        for(int i=1;i<=samples.m;i++){
            typename FILM<T>::SAMPLE& sample=samples(i);
            RENDERING_RAY<T> ray(RAY<VECTOR<T,3> >(camera.position,sample.world_position-camera.position),1,Point_Inside_Object(camera.position));
            dummy_root.recursion_depth=0;dummy_root.current_object=ether;
            sample.radiance=Cast_Ray(ray,dummy_root);
            sample.alpha=ray.ray.semi_infinite?(T)0:(T)1;
            camera.film.Add_Sample(sample);}}
}
//#####################################################################
// Render_Pixel
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Render_Pixel(const VECTOR<int,2>& pixel_index,RAY_TRACER_DEBUG_DATA<T>* debug_data)
{
    RENDERING_RAY<T> dummy_root;
    if(debug_mode&&debug_data){dummy_root.debug_ray=new RENDERING_RAY_DEBUG<T>(dummy_root);debug_data->Set_Ray_Tree(dummy_root.debug_ray);}
    //if(use_adaptive_supersampling) return Compute_Adaptive_Pixel_Color(camera.Pixel_Center(i,j),camera.film->grid.dx,camera.film->grid.dy,0,TV(),TV(),-1,debug_data);
    ARRAY<typename FILM<T>::SAMPLE> samples;
    camera.film.Generate_Stratified_Sample_Points(pixel_index,camera,samples);
    for(int i=1;i<=samples.m;i++){
        typename FILM<T>::SAMPLE& sample=samples(i);
        RENDERING_RAY<T> ray(RAY<VECTOR<T,3> >(camera.position,sample.world_position-camera.position),1,Point_Inside_Object(camera.position));
        dummy_root.recursion_depth=0;dummy_root.current_object=ether;
        sample.radiance=Cast_Ray(ray,dummy_root);
        sample.alpha=ray.ray.semi_infinite?(T)0:(T)1;
        camera.film.Add_Sample(sample);
    }
}
#if 0
//#####################################################################
// Compute_Adaptive_Pixel_Color
//#####################################################################
template<class T> VECTOR<T,3> RENDER_WORLD<T>::
Compute_Adaptive_Pixel_Color(const TV& center,const T dx,const T dy,const int recursion_depth,const TV& precomputed_color_one,const TV& precomputed_color_two,const int precomputed_index,RAY_TRACER_DEBUG_DATA<T>* debug_data)
{
    PHYSBAM_FATAL_ERROR(); // TODO: Fix this


    RENDERING_RAY<T> dummy_root;
    if(debug_mode&&debug_data){dummy_root.debug_ray=new RENDERING_RAY_DEBUG<T>(dummy_root);debug_data->Set_Ray_Tree(dummy_root.debug_ray);}
    
    ARRAY<TV > positions(5);T dx_2=(T)0.5*dx,dy_2=(T)0.5*dy;
    positions(1)=center-dx_2*camera.horizontal_vector-dy_2*camera.vertical_vector;
    positions(2)=center-dx_2*camera.horizontal_vector+dy_2*camera.vertical_vector;
    positions(3)=center+dx_2*camera.horizontal_vector-dy_2*camera.vertical_vector;
    positions(4)=center+dx_2*camera.horizontal_vector+dy_2*camera.vertical_vector;
    positions(5)=center;

    ARRAY<TV > colors(5);
    for(int i=1;i<=5;i++){
        if(i==precomputed_index){colors(i)=precomputed_color_one;continue;}
        if(i==5-precomputed_index){colors(i)=precomputed_color_two;continue;}
        RENDERING_RAY<T> rendering_ray(RAY<VECTOR<T,3> >(camera.position,positions(i)-camera.position),1,Point_Inside_Object(camera.position));
        dummy_root.recursion_depth=0;dummy_root.current_object=ether;
        colors(i)=Cast_Ray(rendering_ray,dummy_root);}
    
    if(recursion_depth==adaptive_supersampling_depth_limit)
        return (T)0.125*colors(1)+(T)0.125*colors(2)+(T)0.125*colors(3)+(T)0.125*colors(4)+(T)0.5*colors(5);

    for(int i=1;i<=4;i++)
        if((colors(i)-colors(5)).Magnitude()<adaptive_supersampling_tolerance) colors(i)=(T)0.5*colors(i)+(T)0.5*colors(5);
        else colors(i)=Compute_Adaptive_Pixel_Color((T)0.5*(positions(i)+positions(5)),dx_2,dy_2,recursion_depth+1,colors(i),colors(5),i,debug_data);

    return (T)0.25*colors(1)+(T)0.25*colors(2)+(T)0.25*colors(3)+(T)0.25*colors(4);
}
#endif
//#####################################################################
// Closest_Intersection
//#####################################################################
template<class T> RENDERING_OBJECT<T>* RENDER_WORLD<T>::
Closest_Intersection(RENDERING_RAY<T>& ray,const int lowest_priority)
{
    const RENDERING_OBJECT<T>* closest_object=0;
    for(int i=1;i<=standard_objects.m;i++)standard_objects(i)->Intersection(ray.ray,lowest_priority,&closest_object);
    return (RENDERING_OBJECT<T>*)closest_object;
}
//#####################################################################
// Point_Inside_Object
//#####################################################################
template<class T> RENDERING_OBJECT<T>* RENDER_WORLD<T>::
Point_Inside_Object(const TV& point)
{
    RENDERING_OBJECT<T> *inside_object=0,*candidate_object=0;
    for(int i=standard_objects.m;i>=1;i--)if(standard_objects(i)->Inside(point,&candidate_object))
        if(!inside_object||candidate_object->priority>=inside_object->priority) inside_object=candidate_object;
    assert(inside_object);return inside_object;
}
//#####################################################################
// Cast_Ray
//#####################################################################
template<class T> VECTOR<T,3> RENDER_WORLD<T>::
Cast_Ray(RENDERING_RAY<T>& ray,const RENDERING_RAY<T>& parent_ray)
{
    if(parent_ray.recursion_depth>=ray_depth_limit-1 || ray.ray_contribution<=ray_contribution_limit) return TV();
    ray.recursion_depth=parent_ray.recursion_depth+1;
    ray.ray_type=RENDERING_RAY<T>::COLOR_RAY;
    RENDERING_OBJECT<T>* object=Closest_Intersection(ray,ray.current_object->priority);
    if(debug_mode && parent_ray.debug_ray){parent_ray.debug_ray->Add_Child(ray);ray.debug_ray->hit_object=(object!=0);}
    if(!object) return Attenuate_Ray(ray,background_shader->Shade_Surface_Using_Direct_Illumination(ray,*ray.current_object,*ray.current_object,*ray.current_object,TV(),TV(1,0,0)));

    TV intersection_point=ray.ray.Point(ray.ray.t_max);
    bool flipped_normal;TV same_side_normal=Same_Side_Normal(ray,intersection_point,*object,flipped_normal);
    if(debug_mode && parent_ray.debug_ray){ray.debug_ray->same_side_normal=same_side_normal;}
    RENDERING_OBJECT<T>* entering_object;
    if(object->support_transparent_overlapping_objects || flipped_normal) entering_object=Point_Inside_Object(intersection_point-same_side_normal*object->small_number*4);
    else entering_object=object;
    TV color=object->material_shader->Shade_Surface_Using_Direct_Illumination(ray,*ray.current_object,*entering_object,*object,intersection_point,same_side_normal);
    if(use_photon_mapping){color+=object->material_shader->Shade_Surface_Using_Indirect_Illumination(ray,*ray.current_object,*entering_object,*object,intersection_point,same_side_normal);}
    return Attenuate_Ray(ray,color);
}
//#####################################################################
// Cast_Ray_For_Photon_Gather
//#####################################################################
template<class T> VECTOR<T,3> RENDER_WORLD<T>::
Cast_Ray_For_Photon_Gather(RENDERING_RAY<T>& ray,const RENDERING_RAY<T>& parent_ray)
{
    if(parent_ray.recursion_depth>=ray_depth_limit-1 || ray.ray_contribution<=ray_contribution_limit) return TV();
    ray.recursion_depth=parent_ray.recursion_depth+1;
    ray.ray_type=RENDERING_RAY<T>::PHOTON_GATHER_RAY;
    RENDERING_OBJECT<T>* object=Closest_Intersection(ray,ray.current_object->priority);
    if(debug_mode && parent_ray.debug_ray){parent_ray.debug_ray->Add_Child(ray);ray.debug_ray->hit_object=(object!=0);}
    if(!object) return Attenuate_Ray(ray,background_shader->Shade_Surface_Using_Direct_Illumination(ray,*ray.current_object,*ray.current_object,*ray.current_object,TV(),TV(1,0,0)));

    TV intersection_point=ray.ray.Point(ray.ray.t_max);
    bool flipped_normal;TV same_side_normal=Same_Side_Normal(ray,intersection_point,*object,flipped_normal);
    if(debug_mode && parent_ray.debug_ray){ray.debug_ray->same_side_normal=same_side_normal;}
    RENDERING_OBJECT<T>* entering_object;
    if(object->support_transparent_overlapping_objects || flipped_normal) entering_object=Point_Inside_Object(intersection_point-same_side_normal*object->small_number*4);
    else entering_object=object;
    TV color=object->material_shader->Shade_Surface_Using_Approximate_Full_Illumination(ray,*ray.current_object,*entering_object,*object,intersection_point,same_side_normal);
    return Attenuate_Ray(ray,color);
}
//#####################################################################
// Incident_Light 
//#####################################################################
template<class T> VECTOR<T,3> RENDER_WORLD<T>::
Incident_Light(RENDERING_RAY<T> ray,const RENDERING_LIGHT<T>& light,const RENDERING_RAY<T>& full_ray,const RENDERING_RAY<T>& parent_ray)
{
    if(parent_ray.recursion_depth>=ray_depth_limit-1 || ray.ray_contribution<=ray_contribution_limit) return TV();
    RENDERING_OBJECT<T>* object=0;
    if(light.casts_shadows){
        ray.recursion_depth=parent_ray.recursion_depth+1;
        ray.ray_type=RENDERING_RAY<T>::SHADOW_RAY;
        object=Closest_Intersection(ray,ray.current_object->priority);
        if(debug_mode && parent_ray.debug_ray){parent_ray.debug_ray->Add_Child(ray);ray.debug_ray->hit_object=(object!=0);}}
    TV incident_light_color;
    if(object){
        TV intersection_point=ray.ray.Point(ray.ray.t_max);
        bool flipped_normal;TV same_side_normal=Same_Side_Normal(ray,intersection_point,*object,flipped_normal);
        if(debug_mode && parent_ray.debug_ray){ray.debug_ray->same_side_normal=same_side_normal;}
        RENDERING_OBJECT<T>* entering_object;
        if(object->support_transparent_overlapping_objects||flipped_normal) entering_object=Point_Inside_Object(intersection_point-same_side_normal*object->small_number*4);
        else entering_object=object;
        incident_light_color=object->material_shader->Shade_Light_Ray(ray,*ray.current_object,*entering_object,*object,intersection_point,same_side_normal,light,full_ray);}
    else incident_light_color=light.Emitted_Light(full_ray);
    return Attenuate_Light_Ray(ray,light,incident_light_color);
}
//#####################################################################
// Incident_Light 
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Cast_Photon(RENDERING_RAY<T>& incoming_ray,const RENDERING_RAY<T>& parent_ray,const TV& power,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,
            const int diffuse_bounces,const int specular_bounces)
{
    if(parent_ray.recursion_depth>=photon_depth_limit-1)return;
    incoming_ray.recursion_depth=parent_ray.recursion_depth+1;
    incoming_ray.ray_type=RENDERING_RAY<T>::PHOTON_RAY;
    RENDERING_OBJECT<T>* scattering_medium_object;RENDERING_OBJECT<T>* object;
    RENDERING_RAY<T> working_ray=incoming_ray;
    TV current_power=power;
    T fixed_step_size=(T).01;

    while(true){
        object=Closest_Intersection(working_ray,working_ray.current_object->priority);
        // if object is less than "fixed step size away", break
        if(object && (working_ray.ray.Point(working_ray.ray.t_max)-working_ray.ray.endpoint).Magnitude()<=fixed_step_size) break;

        scattering_medium_object=0;
        for(int i=1;i<=volumetric_objects.m;i++)
            if(volumetric_objects(i)->volumetric_shader->supports_photon_mapping&&volumetric_objects(i)->Inside(working_ray.ray.endpoint))
                scattering_medium_object=volumetric_objects(i);
        if(scattering_medium_object){
            // passing through a scattering medium. handle event. if absorb, stop.  otherwise, go back to beginning of loop.
            bool continue_scattering=scattering_medium_object->volumetric_shader->Scatter_Photon_Ray(*scattering_medium_object,working_ray,current_power,type,fixed_step_size);
            if(!continue_scattering)return;}
        else{
            // not initially in a scattering medium, test if we hit it before object intersection
            for(int i=1;i<=volumetric_objects.m;i++)
                if(volumetric_objects(i)->volumetric_shader->supports_photon_mapping&&volumetric_objects(i)->Intersection(working_ray.ray))
                    scattering_medium_object=volumetric_objects(i);
            if(scattering_medium_object){
                // passes through a scattering medium, just not right away, so change ray endpoint to that point
                TV intersection_point=working_ray.ray.Point(working_ray.ray.t_max);
                TV normal=scattering_medium_object->Normal(intersection_point,working_ray.ray.aggregate_id);
                if(TV::Dot_Product(normal,working_ray.ray.direction)>0)normal*=-1;
                working_ray.ray.endpoint=intersection_point-normal*scattering_medium_object->small_number;
                working_ray.ray.t_max=0;working_ray.ray.semi_infinite=true;}
            else break;}}

    if(object){
        TV intersection_point=working_ray.ray.Point(working_ray.ray.t_max);incoming_ray.ray.aggregate_id=working_ray.ray.aggregate_id;
        bool flipped_normal;TV same_side_normal=Same_Side_Normal(working_ray,intersection_point,*object,flipped_normal);
        RENDERING_OBJECT<T>* entering_object;
        if(object->support_transparent_overlapping_objects||flipped_normal) entering_object=Point_Inside_Object(intersection_point-same_side_normal*object->small_number*2);
        else entering_object=object;
        bool should_throw=false;
        TV attenuated_power=Attenuate_Photon(working_ray,current_power,should_throw);
        if(!should_throw)
            object->material_shader->Receive_Photon(incoming_ray,*incoming_ray.current_object,*entering_object,*object,intersection_point,same_side_normal,
                attenuated_power,type,diffuse_bounces,specular_bounces);}
}
//#####################################################################
// Attenuate_Ray
//#####################################################################
template<class T> VECTOR<T,3> RENDER_WORLD<T>::
Attenuate_Ray(RENDERING_RAY<T>& ray,const TV& color)
{
    TV attenuated_color=color;
    for(int i=1;i<=volumetric_objects.m;i++)
        if(volumetric_objects(i)->priority>=ray.current_object->priority)
            attenuated_color=volumetric_objects(i)->volumetric_shader->Attenuate_Color(ray,*volumetric_objects(i),attenuated_color);
    return attenuated_color;
}
//#####################################################################
// Attenuate_Light_Ray
//#####################################################################
template<class T> VECTOR<T,3> RENDER_WORLD<T>::
Attenuate_Light_Ray(RENDERING_RAY<T>& ray,const RENDERING_LIGHT<T>& light,const TV& color)
{
    TV attenuated_light_color=color;
    for(int i=1;i<=volumetric_objects.m;i++) 
        if(volumetric_objects(i)->priority>=ray.current_object->priority)
            attenuated_light_color=volumetric_objects(i)->volumetric_shader->Attenuate_Light(ray,*volumetric_objects(i),light,attenuated_light_color);
    return attenuated_light_color;
}
//#####################################################################
// Attenuate_Photon
//#####################################################################
template<class T> VECTOR<T,3> RENDER_WORLD<T>::
Attenuate_Photon(RENDERING_RAY<T>& ray,const TV& photon_power,bool& should_throw)
{
    TV attenuated_photon_power=photon_power;
    should_throw=false;
    for(int i=1;i<=volumetric_objects.m;i++)
        if(volumetric_objects(i)->priority>=ray.current_object->priority){
            attenuated_photon_power=volumetric_objects(i)->volumetric_shader->Attenuate_Photon(ray,*volumetric_objects(i),attenuated_photon_power,should_throw);
            if(should_throw) break;
        }
    return attenuated_photon_power;
}
//#####################################################################
// Prepare_For_Forward_Ray_Trace
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Prepare_For_Forward_Ray_Trace()
{
    for(int i=1;i<=standard_objects.m;i++)standard_objects(i)->Preprocess_Efficiency_Structures(*this);
    for(int i=1;i<=volumetric_objects.m;i++)volumetric_objects(i)->Preprocess_Efficiency_Structures(*this);

    if(use_photon_mapping){ // TODO: adjust photons per light based on individual light power
        // count photon mapped lights
        int global_photon_emitting_lights=0,caustic_photon_emitting_lights=0,volume_photon_emitting_lights=0;
        for(int i=1;i<=lights.m;i++){if(lights(i)->supports_global_photon_mapping)global_photon_emitting_lights++;if(lights(i)->supports_caustic_photon_mapping)caustic_photon_emitting_lights++;
            if(lights(i)->supports_volume_photon_mapping)volume_photon_emitting_lights++;}
        int global_photons_per_light=global_photon_emitting_lights?global_photon_map.Max_Number_Of_Photons()/global_photon_emitting_lights:0;
        int global_photons_remainder=global_photon_emitting_lights?global_photon_map.Max_Number_Of_Photons()%global_photon_emitting_lights:0;
        int caustic_photons_per_light=caustic_photon_emitting_lights?caustic_photon_map.Max_Number_Of_Photons()/caustic_photon_emitting_lights:0;
        int caustic_photons_remainder=caustic_photon_emitting_lights?caustic_photon_map.Max_Number_Of_Photons()%caustic_photon_emitting_lights:0;
        int volume_photons_per_light=volume_photon_emitting_lights?volume_photon_map.Max_Number_Of_Photons()/volume_photon_emitting_lights:0;
        int volume_photons_remainder=volume_photon_emitting_lights?volume_photon_map.Max_Number_Of_Photons()%volume_photon_emitting_lights:0;

        // shoot the photons
        {LOG::SCOPE scope("shooting photons","Shooting Photons...%d global %d caustic %d volume\n",global_photon_map.Max_Number_Of_Photons(),caustic_photon_map.Max_Number_Of_Photons(),
            volume_photon_map.Max_Number_Of_Photons());
        RENDERING_RAY<T> dummy_ray;dummy_ray.current_object=ether;dummy_ray.recursion_depth=0;
        for(int i=1;i<=lights.m;i++){
            LOG::SCOPE scope("light","light %d",i);
            if(lights(i)->supports_global_photon_mapping){
                LOG::Time("global photons");
                int global_photon_budget=global_photons_per_light+((i==lights.m)?global_photons_remainder:0);
                Shoot_Photons_From_Light_For_Specific_Map(dummy_ray,i,global_photon_map,global_photon_budget,PHOTON_MAP<T>::GLOBAL_PHOTON_MAP);}
            if(lights(i)->supports_caustic_photon_mapping){
                LOG::Time("caustic photons");
                int caustic_photon_budget=caustic_photons_per_light+((i==lights.m)?caustic_photons_remainder:0);
                Shoot_Photons_From_Light_For_Specific_Map(dummy_ray,i,caustic_photon_map,caustic_photon_budget,PHOTON_MAP<T>::CAUSTIC_PHOTON_MAP);}
            if(lights(i)->supports_volume_photon_mapping){
                LOG::Time("volume photons");
                int volume_photon_budget=volume_photons_per_light+((i==lights.m)?volume_photons_remainder:0);
                Shoot_Photons_From_Light_For_Specific_Map(dummy_ray,i,volume_photon_map,volume_photon_budget,PHOTON_MAP<T>::VOLUME_PHOTON_MAP);}}}

        // construct balanced KD-Tree's
        LOG::Time("balancing kd-trees");
        global_photon_map.Prepare_Photon_Map_For_Rendering();
        caustic_photon_map.Prepare_Photon_Map_For_Rendering();
        volume_photon_map.Prepare_Photon_Map_For_Rendering();
        LOG::Stop_Time();
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        if(use_irradiance_cache)irradiance_cache.Enable_Cache(global_photon_map.bounding_box,irradiance_cache_max_distance);
#endif
    }
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        Prepare_BSSRDF_Octree();
#endif
}
//#####################################################################
// Shoot_Photons_From_Light_For_Specific_Map
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Shoot_Photons_From_Light_For_Specific_Map(RENDERING_RAY<T>& parent_ray,int light_index,PHOTON_MAP<T>& map,int photon_budget,typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type)
{
    map.Begin_Light_Emission(photon_budget);
    int number_of_photons_emitted=lights(light_index)->Emit_Photons(parent_ray,map,type);
    map.End_Light_Emission(number_of_photons_emitted);
}
//#####################################################################
// Use_Photon_Mapping
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Use_Photon_Mapping(const int global_photons,const int caustic_photons,const int volume_photons,const T max_photon_distance,const int number_of_photons_for_estimate)
{
    use_photon_mapping=true;
    global_photon_map.Resize_Photon_Map(global_photons);
    caustic_photon_map.Resize_Photon_Map(caustic_photons);
    volume_photon_map.Resize_Photon_Map(volume_photons);
    this->max_photon_distance=max_photon_distance;
    this->number_of_photons_for_estimate=number_of_photons_for_estimate;
}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Use_Irradiance_Cache
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Use_Irradiance_Cache(const bool use_irradiance_cache_input,const T irradiance_cache_max_distance_input,const int irradiance_cache_samples_input)
{
    irradiance_cache.Set_Final_Gather_Samples(irradiance_cache_samples_input); // irradiance cache uses this for final gather
    use_irradiance_cache=use_irradiance_cache_input;irradiance_cache_max_distance=irradiance_cache_max_distance_input;
}
#endif
//#####################################################################
// Use_Adaptive_Sampling
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Use_Adaptive_Supersampling(const bool use_adaptive_supersampling_input,const T adaptive_supersampling_tolerance_input,const int adaptive_supersampling_depth_limit_input)
{
    use_adaptive_supersampling=use_adaptive_supersampling_input;adaptive_supersampling_tolerance=adaptive_supersampling_tolerance_input;
    adaptive_supersampling_depth_limit=adaptive_supersampling_depth_limit_input;
}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Prepare_BSSRDF_Octree
//#####################################################################
template<class T> void RENDER_WORLD<T>::
Prepare_BSSRDF_Octree()
{
    for(int i=1;i<=standard_objects.m;i++) if(standard_objects(i)->bssrdf_shader) standard_objects(i)->Generate_BSSRDF_Tree(*this);
}
#endif
//#####################################################################
template class RENDER_WORLD<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDER_WORLD<double>;
#endif
