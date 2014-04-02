//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Tessellation/CYLINDER_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const CYLINDER<T>& cylinder,const int resolution_height,const int resolution_radius)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    surface->Initialize_Cylinder_Mesh_And_Particles(resolution_height,resolution_radius,cylinder.height,cylinder.radius,true);
    FRAME<TV> frame(cylinder.plane1.x1,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),cylinder.plane2.x1-cylinder.plane1.x1));
    particles.X=frame*particles.X;
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const CYLINDER<float>&,const int,const int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const CYLINDER<double>&,const int,const int);
#endif
}
}
