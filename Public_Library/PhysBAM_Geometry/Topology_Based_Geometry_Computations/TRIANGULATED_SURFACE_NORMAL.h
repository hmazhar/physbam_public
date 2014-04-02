//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_NORMAL__
#define __TRIANGULATED_SURFACE_NORMAL__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T> VECTOR<T,3>
Normal(const TRIANGULATED_SURFACE<T>& ts,const VECTOR<T,3>& location,const int aggregate)
{  
    typedef VECTOR<T,3> TV;
    assert(aggregate >= 1 && aggregate <= ts.mesh.elements.m);

    if(ts.use_vertex_normals){
        TV normal1,normal2,normal3;
        int node1,node2,node3;ts.mesh.elements(aggregate).Get(node1,node2,node3);
        if(ts.avoid_normal_interpolation_across_sharp_edges) (*ts.face_vertex_normals)(aggregate).Get(normal1,normal2,normal3);
        else{normal1=(*ts.vertex_normals)(node1);normal2=(*ts.vertex_normals)(node2);normal3=(*ts.vertex_normals)(node3);}
        TV weights=TRIANGLE_3D<T>::Barycentric_Coordinates(location,ts.particles.X(node1),ts.particles.X(node2),ts.particles.X(node3));
        return (weights.x*normal1+weights.y*normal2+weights.z*normal3).Normalized();}
    else return ts.Face_Normal(aggregate);
}
}
}
#endif
