//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MESH_OBJECT_PRUNE__
#define __MESH_OBJECT_PRUNE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>

namespace PhysBAM{
template<class TV,class T_MESH> class MESH_OBJECT;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Discard_Valence_Zero_Particles_And_Renumber
//#####################################################################
template<class TV,class T_MESH>
void Discard_Valence_Zero_Particles_And_Renumber(MESH_OBJECT<TV,T_MESH>& mo,ARRAY<int>& condensation_mapping)
{
    // mark which nodes are used
    ARRAY<bool> node_is_used(mo.mesh.number_nodes);
    for(int t=1;t<=mo.mesh.elements.m;t++){
        INDIRECT_ARRAY<ARRAY<bool>,typename T_MESH::ELEMENT_TYPE&> subset=node_is_used.Subset(mo.mesh.elements(t));
        ARRAYS_COMPUTATIONS::Fill(subset,true);}
    
    // make condensation mapping
    condensation_mapping.Resize(mo.mesh.number_nodes,false,false);ARRAYS_COMPUTATIONS::Fill(condensation_mapping,0);
    int counter=0;
    for(int t=1;t<=mo.mesh.number_nodes;t++) if(node_is_used(t)) condensation_mapping(t)=++counter;
    
    // make new triangle mesh
    mo.mesh.number_nodes=counter;
    for(int t=1;t<=mo.mesh.elements.m;t++)
        mo.mesh.elements(t)=condensation_mapping.Subset(mo.mesh.elements(t));
    
    // do particles same way
    for(int p=1;p<=condensation_mapping.m;p++) if(!condensation_mapping(p)) mo.particles.array_collection->Add_To_Deletion_List(p);
    for(int p=condensation_mapping.m+1;p<=mo.particles.array_collection->Size();p++) mo.particles.array_collection->Add_To_Deletion_List(p);
    mo.particles.array_collection->Delete_Elements_On_Deletion_List(true);mo.particles.array_collection->Compact();

    mo.Refresh_Auxiliary_Structures();
}

}
}
#endif
