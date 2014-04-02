//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MESH_OBJECT_UNION__
#define __MESH_OBJECT_UNION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>

namespace PhysBAM{
namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS{
//#####################################################################
// Function Union_Mesh_Objects_Relatively
//#####################################################################
template<class TV,class T_OBJECT>
void Union_Mesh_Objects_Relatively(T_OBJECT *object,const ARRAY<T_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames)
{
    GEOMETRY_PARTICLES<TV>& particles=object->particles;
    particles.array_collection->Clean_Memory();
    object->mesh.elements.Remove_All();
    // resize
    {int total_particles=0,total_elements=0;
    for(int i=1;i<=object_list.m;i++){
        total_particles+=object_list(i)->particles.array_collection->Size();
        total_elements+=object_list(i)->mesh.elements.m;}
    particles.array_collection->Preallocate(total_particles);object->mesh.elements.Preallocate(total_elements);}
    // copy
    for(int i=1;i<=object_list.m;i++){
        int particle_offset=particles.array_collection->Size();
        particles.array_collection->Add_Arrays(*object_list(i)->particles.array_collection);
        particles.array_collection->Append(*object_list(i)->particles.array_collection);
        for(int p=1;p<=object_list(i)->particles.array_collection->Size();p++){int p2=p+particle_offset;
            particles.X(p2)=relative_frames(i)*particles.X(p2);
            if(particles.store_velocity) particles.V(p2)=relative_frames(i).r.Rotate(particles.V(p2));}
        object->mesh.elements.Append_Elements(object_list(i)->mesh.elements+particle_offset);}
    object->Update_Number_Nodes();
}

}
}
#endif
