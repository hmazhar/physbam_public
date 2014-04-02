//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Tessellation/ORIENTED_BOX_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const ORIENTED_BOX<VECTOR<T,3> >& box)
{
    typedef VECTOR<T,3> TV;
    MATRIX<T,3> directions=box.edges;TV lengths;
    for(int i=1;i<=3;++i) lengths[i]=directions.Column(i).Normalize();
    BOX<TV> aligned_box(TV(),lengths);
    TRIANGULATED_SURFACE<T>* surface=Generate_Triangles(aligned_box);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    FRAME<TV> frame(box.corner,ROTATION<TV>(directions));
    particles.X=frame*particles.X;
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const ORIENTED_BOX<VECTOR<float,3> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const ORIENTED_BOX<VECTOR<double,3> >&);
#endif
}
}
