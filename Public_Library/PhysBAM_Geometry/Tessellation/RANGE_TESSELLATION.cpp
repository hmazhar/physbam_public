//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RANGE<VECTOR<T,3> >& box)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(8);
    particles.X(1)=TV(box.min_corner.x,box.min_corner.y,box.max_corner.z);particles.X(2)=TV(box.max_corner.x,box.min_corner.y,box.max_corner.z);
    particles.X(3)=TV(box.max_corner.x,box.max_corner.y,box.max_corner.z);particles.X(4)=TV(box.min_corner.x,box.max_corner.y,box.max_corner.z);
    particles.X(5)=TV(box.min_corner.x,box.min_corner.y,box.min_corner.z);particles.X(6)=TV(box.max_corner.x,box.min_corner.y,box.min_corner.z);
    particles.X(7)=TV(box.max_corner.x,box.max_corner.y,box.min_corner.z);particles.X(8)=TV(box.min_corner.x,box.max_corner.y,box.min_corner.z);
    TRIANGLE_MESH& mesh=surface->mesh;mesh.number_nodes=8;mesh.elements.Preallocate(12);
    mesh.elements.Append(VECTOR<int,3>(1,4,8));mesh.elements.Append(VECTOR<int,3>(8,5,1));
    mesh.elements.Append(VECTOR<int,3>(2,6,7));mesh.elements.Append(VECTOR<int,3>(7,3,2));
    mesh.elements.Append(VECTOR<int,3>(5,2,1));mesh.elements.Append(VECTOR<int,3>(2,5,6));
    mesh.elements.Append(VECTOR<int,3>(7,8,3));mesh.elements.Append(VECTOR<int,3>(4,3,8));
    mesh.elements.Append(VECTOR<int,3>(5,8,7));mesh.elements.Append(VECTOR<int,3>(7,6,5));
    mesh.elements.Append(VECTOR<int,3>(1,2,3));mesh.elements.Append(VECTOR<int,3>(1,3,4));
    return surface;
}
//#####################################################################
}
template TRIANGULATED_SURFACE<float>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<float,3> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<double,3> >&);
#endif
}
