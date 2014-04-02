//#####################################################################
// Copyright 2004, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_KD_TREE_ACCELERATOR
//#####################################################################
#ifndef __RENDERING_KD_TREE_ACCELERATOR__
#define __RENDERING_KD_TREE_ACCELERATOR__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_KD_TREE_ACCELERATOR:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::material_shader;using RENDERING_OBJECT<T>::volumetric_shader;

    ARRAY<RENDERING_OBJECT<T>*> objects;
    //KD_TREE_3D<T,RENDERING_OBJECT<T>*> kd_tree;
    ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> > primitives;
    RANGE<TV> bounding_box;

    RENDERING_KD_TREE_ACCELERATOR()
    {}

    void Add_Object(RENDERING_OBJECT<T>* object)
    {// so this get's into the right bin in render world
    if(object->material_shader)material_shader=object->material_shader; 
    if(object->volumetric_shader)volumetric_shader=object->volumetric_shader;
    objects->Append(object);}

    bool Intersection(RAY<TV>& ray,const int lowest_priority,RENDERING_OBJECT<T>** intersected_object)const PHYSBAM_OVERRIDE
    {bool hit=false;
    for(int i=1;i<=primitives.m;i++){
        if(primitives(i).object->priority>=lowest_priority){
            bool primitive_i_hit=primitives(i).object->Intersection(ray,primitives(i).aggregate_id) PHYSBAM_OVERRIDE;
            if(primitive_i_hit){hit=true;*intersected_object=primitives(i).object;}}}
    return hit;}

    void Preprocess_Efficiency_Structures() PHYSBAM_OVERRIDE
    {for(int i=1;i<=objects.m;i++)objects(i)->Get_Aggregate_World_Space_Bounding_Boxes(primitives);
    if(primitives.m>0){
        // construct total bounding box...
        bounding_box.Reset_Bounds(primitives(1).bounding_box);
        for(int i=2;i<=objects.m;i++)bounding_box.Enlarge_To_Include_Box(primitives(i).bounding_box);}}

    bool Inside(const TV& location,RENDERING_OBJECT<T>** intersected_object) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=objects.m;i++){
        if(objects(i)->support_transparent_overlapping_objects&&objects(i)->Inside(location)){
            *intersected_object=(RENDERING_OBJECT<T>*)this;return true;}}
    return false;}

//#####################################################################
};   
}
#endif

