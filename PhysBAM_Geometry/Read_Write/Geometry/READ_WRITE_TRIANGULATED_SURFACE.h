//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TRIANGULATED_SURFACE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TRIANGULATED_SURFACE__
#define __READ_WRITE_TRIANGULATED_SURFACE__

#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

void Register_Read_Write_Triangulated_Surface();

template<class RW,class T>
class Read_Write<TRIANGULATED_SURFACE<T>,RW>:public Read_Write<MESH_OBJECT<VECTOR<T,3>,TRIANGLE_MESH>,RW>
{
public:
    static void Print_Statistics(std::ostream& output,TRIANGULATED_SURFACE<T>& object,const T thickness_over_2=(T)1e-8)
    {if(object.mesh.number_nodes!=object.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    int index;object.Update_Bounding_Box();
    if(!object.mesh.incident_elements) object.mesh.Initialize_Incident_Elements();

    output<<"triangles = "<<object.mesh.elements.m<<std::endl;
    output<<"particles = "<<object.particles.array_collection->Size()<<std::endl;
    {int particles_touched=0;for(int p=1;p<=object.particles.array_collection->Size();p++) if((*object.mesh.incident_elements)(p).m) particles_touched++;
    output<<"particles touched = "<<particles_touched<<std::endl;}
    output<<"bounding box = "<<*object.bounding_box<<std::endl;
    if(object.particles.store_velocity){
        int index=ARRAYS_COMPUTATIONS::Arg_Maximum_Magnitude(object.particles.V);
        output<<"max_speed = "<<object.particles.V(index).Magnitude()<<" ("<<index<<")"<<std::endl;}
    output<<"total area = "<<object.Total_Area()<<std::endl;
    output<<"min area = "<<object.Minimum_Area(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"min altitude = "<<object.Minimum_Altitude(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"min edge length = "<<object.Minimum_Edge_Length(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max edge length = "<<object.Maximum_Edge_Length(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"min angle = "<<object.Minimum_Angle(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max angle = "<<object.Maximum_Angle(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"ave min angle = "<<object.Average_Minimum_Angle()<<std::endl;
    output<<"ave max angle = "<<object.Average_Maximum_Angle()<<std::endl;
    output<<"max aspect ratio = "<<object.Maximum_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"ave aspect ratio = "<<object.Average_Aspect_Ratio()<<std::endl;
    if(object.Check_For_Self_Intersection(thickness_over_2)) output<<"found self intersections"<<std::endl;else output<<"no self intersections"<<std::endl;}
};
}
#endif
#endif
