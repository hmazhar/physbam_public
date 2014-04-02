//#####################################################################
// Copyright 2004, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_GRID_ACCELERATOR
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_GRID_ACCELERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Class INTERSECTION_MAP_HELPER
//#####################################################################
template<class T> class INTERSECTION_MAP_HELPER {
public:
    const RENDERING_UNIFORM_GRID_ACCELERATOR<T>* rendering_uniform_accelerator;
    int lowest_priority;
    const RENDERING_OBJECT<T>** intersected_object;RAY<VECTOR<T,3> > working_ray;
    
    bool Callback(RAY<VECTOR<T,3> >& ray,const ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>*>& primitives,const T cell_t_max)
    {const RENDERING_OBJECT<T>* closest_object=0;
    for(int primitive_id=1;primitive_id<=primitives.m;primitive_id++){
        RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>& primitive=*primitives(primitive_id);
        if(primitive.operation==rendering_uniform_accelerator->operation){ // already have done this intersection
            if(primitive.hit_aggregate_id!=-1&&primitive.hit_t<=cell_t_max&&primitive.hit_t<=working_ray.t_max){
                closest_object=primitive.object;working_ray.t_max=primitive.hit_t;working_ray.aggregate_id=primitive.hit_aggregate_id;}}
        else if(primitive.object->Intersection(working_ray,primitive.aggregate_id)){ // got intersection
              if(working_ray.t_max<=cell_t_max)closest_object=primitive.object; // intersection was in cell
              else{ // intersection was not in current cell, save for later
                  primitive.hit_t=working_ray.t_max;primitive.hit_aggregate_id=working_ray.aggregate_id;primitive.operation=rendering_uniform_accelerator->operation;}}
        else{ // no intersection, save for later
            primitive.hit_aggregate_id=-1;primitive.operation=rendering_uniform_accelerator->operation;}}
    if(closest_object){*intersected_object=closest_object;ray.Restore_Intersection_Information(working_ray);return true;}
    else return false;}
};
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_UNIFORM_GRID_ACCELERATOR<T>::
Intersection(RAY<VECTOR<T,3> >& ray,const int lowest_priority,const RENDERING_OBJECT<T>** intersected_object) const
{
    INTERSECTION_MAP_HELPER<T> helper;helper.rendering_uniform_accelerator=this;helper.lowest_priority=lowest_priority;helper.intersected_object=intersected_object;
    helper.working_ray=RAY<VECTOR<T,3> >(ray);
    bool intersect_value=uniform_grid.template Map_Intersection<INTERSECTION_MAP_HELPER<T>*>(ray,&helper);
    operation++; // throw away all saved mailboxes
    return intersect_value;
}
//#####################################################################
template class RENDERING_UNIFORM_GRID_ACCELERATOR<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDERING_UNIFORM_GRID_ACCELERATOR<double>;
#endif
