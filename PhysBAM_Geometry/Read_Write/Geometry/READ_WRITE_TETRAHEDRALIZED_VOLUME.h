//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TETRAHEDRALIZED_VOLUME
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TETRAHEDRALIZED_VOLUME__
#define __READ_WRITE_TETRAHEDRALIZED_VOLUME__

#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{

void Register_Read_Write_Tetrahedralized_Volume();

template<class RW,class T>
class Read_Write<TETRAHEDRALIZED_VOLUME<T>,RW>:public Read_Write<MESH_OBJECT<VECTOR<T,3>,TETRAHEDRON_MESH>,RW>
{
public:
    static void Print_Statistics(std::ostream& output,TETRAHEDRALIZED_VOLUME<T>& object)
    {if(object.mesh.number_nodes!=object.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    object.Update_Bounding_Box();
    if(!object.mesh.segment_mesh) object.mesh.Initialize_Segment_Mesh();
    if(!object.mesh.incident_elements) object.mesh.Initialize_Incident_Elements();
    if(!object.mesh.adjacent_elements) object.mesh.Initialize_Adjacent_Elements();
    if(!object.mesh.boundary_mesh) object.mesh.Initialize_Boundary_Mesh();

    output<<"tetrahedrons = "<<object.mesh.elements.m<<std::endl;
    output<<"particles = "<<object.particles.array_collection->Size()<<std::endl;
    {int particles_touched=0;for(int p=1;p<=object.particles.array_collection->Size();p++) if((*object.mesh.incident_elements)(p).m) particles_touched++;
    output<<"particles touched = "<<particles_touched<<std::endl;}
    output<<"bounding box = "<<*object.bounding_box<<std::endl;
    if(object.particles.store_velocity){
        int index=ARRAYS_COMPUTATIONS::Arg_Maximum_Magnitude(object.particles.V);
        output<<"max_speed = "<<object.particles.V(index).Magnitude()<<" ("<<index<<")"<<std::endl;}
    int index;
    output<<"total volume = "<<object.Total_Volume()<<std::endl;
    output<<"max_aspect_ratio = "<<object.Maximum_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max_boundary_aspect_ratio = "<<object.Maximum_Boundary_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max_interior_aspect_ratio = "<<object.Maximum_Interior_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"avg_boundary_aspect_ratio = "<<object.Average_Boundary_Aspect_Ratio()<<std::endl;
    output<<"avg_interior_aspect_ratio = "<<object.Average_Interior_Aspect_Ratio()<<std::endl;
    output<<"min_volume = "<<object.Minimum_Volume(&index)<<std::endl;
    output<<"min_angle = "<<180/pi*object.Minimum_Angle()<<std::endl;
    output<<"max_angle = "<<180/pi*object.Maximum_Angle()<<std::endl;
    output<<"min_dihedral_angle = "<<180/pi*object.Minimum_Dihedral_Angle()<<std::endl;
    output<<"max_dihedral_angle = "<<180/pi*object.Maximum_Dihedral_Angle()<<std::endl;
    output<<"min_edge_length = "<<object.Minimum_Edge_Length()<<std::endl;
    output<<"max_edge_length = "<<object.Maximum_Edge_Length()<<std::endl;
    output<<"min_altitude = "<<object.Minimum_Altitude()<<std::endl;

    ARRAY<int> nonmanifold_nodes;object.mesh.boundary_mesh->Non_Manifold_Nodes(nonmanifold_nodes);
    output<<nonmanifold_nodes.m<<" nonmanifold nodes = "<<nonmanifold_nodes;}

};
}
#endif
#endif
