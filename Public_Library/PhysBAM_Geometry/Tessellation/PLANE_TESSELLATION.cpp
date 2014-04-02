//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLANE_TESSELLATION
//##################################################################### 
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Tessellation/PLANE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const PLANE<T>& plane)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    surface->particles.array_collection->Add_Elements(4);
    TV u=(T)1000*plane.normal.Unit_Orthogonal_Vector(),v=TV::Cross_Product(plane.normal,u);
    surface->particles.X(1)=plane.x1+u+v;
    surface->particles.X(2)=plane.x1-u+v;
    surface->particles.X(3)=plane.x1-u-v;
    surface->particles.X(4)=plane.x1+u-v;
    surface->mesh.elements.Append(VECTOR<int,3>(1,2,3));
    surface->mesh.elements.Append(VECTOR<int,3>(1,3,4));
    surface->Update_Number_Nodes();
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const PLANE<float>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const PLANE<double>&);
#endif
}
}
