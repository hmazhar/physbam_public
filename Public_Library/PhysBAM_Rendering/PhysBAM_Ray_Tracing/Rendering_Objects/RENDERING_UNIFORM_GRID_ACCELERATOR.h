//#####################################################################
// Copyright 2004, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_GRID_ACCELERATOR
//#####################################################################
#ifndef __RENDERING_UNIFORM_GRID_ACCELERATOR__
#define __RENDERING_UNIFORM_GRID_ACCELERATOR__

#include <PhysBAM_Geometry/Spatial_Acceleration/UNIFORM_BOX_PARTITION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_UNIFORM_GRID_ACCELERATOR:public RENDERING_OBJECT<T>
{
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
public:
    using RENDERING_OBJECT<T>::material_shader;using RENDERING_OBJECT<T>::volumetric_shader;

    ARRAY<RENDERING_OBJECT<T>*> objects;
    UNIFORM_BOX_PARTITION<T,RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>*> uniform_grid;
    ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> > primitives;
    mutable int operation;
    
    RENDERING_UNIFORM_GRID_ACCELERATOR()
        :operation(1)
    {}

    void Add_Object(RENDERING_OBJECT<T>* object) // so this get's into the right bin in render world
    {if(object->material_shader)material_shader=object->material_shader; 
    if(object->volumetric_shader)volumetric_shader=object->volumetric_shader;
    objects.Append(object);}
        
    void Preprocess_Efficiency_Structures(RENDER_WORLD<T>& world) PHYSBAM_OVERRIDE
    {for(int i=1;i<=objects.m;i++){objects(i)->Get_Aggregate_World_Space_Bounding_Boxes(primitives);objects(i)->Preprocess_Efficiency_Structures(world);}
    ARRAY<PAIR<RANGE<VECTOR<T,3> >,RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>*> > box_input;
    for(int i=1;i<=primitives.m;i++) box_input.Append(PAIR<RANGE<VECTOR<T,3> >,RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>*>(primitives(i).world_bounding_box,&primitives(i)));
    uniform_grid.Initialize(box_input);}
    
    bool Inside(const VECTOR<T,3>& location,RENDERING_OBJECT<T>** intersected_object) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=objects.m;i++) if(objects(i)->support_transparent_overlapping_objects&&objects(i)->Inside(location)){*intersected_object=(RENDERING_OBJECT<T>*)this;return true;}
    return false;}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {return 0;};

//#####################################################################
    bool Intersection(RAY<VECTOR<T,3> >& ray,const int lowest_priority,const RENDERING_OBJECT<T>** intersected_object) const;
//#####################################################################
};   
}
#endif

