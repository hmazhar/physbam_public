//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGLE_SUBDIVISION_ROOT_THREE__
#define __TRIANGLE_SUBDIVISION_ROOT_THREE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Apply_Root_Three_Subdivision
//#####################################################################
template<class TV> void
Apply_Root_Three_Subdivision(TRIANGLE_SUBDIVISION& ts,ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==ts.start_index_for_new_nodes-1+ts.triangle_mesh.elements.m);
    if(!ts.triangle_mesh.neighbor_nodes){ts.delete_neighbor_nodes=true;ts.triangle_mesh.Initialize_Neighbor_Nodes();}
    for(int p=1;p<=ts.triangle_mesh.number_nodes;p++){
        int n=(*ts.triangle_mesh.neighbor_nodes)(p).m;if(n<3)continue;T one_over_n=(T)1/n;
        T alpha=T(2./9)*(2-cos(T(2*pi)*one_over_n));
        TV sum=base_values((*ts.triangle_mesh.neighbor_nodes)(p)(1));for(int i=2;i<=n;i++)sum+=base_values((*ts.triangle_mesh.neighbor_nodes)(p)(i));
        subdivided_values(p)=(1-alpha)*base_values(p)+alpha*one_over_n*sum;}
    for(int t=1;t<=ts.triangle_mesh.elements.m;t++){
        int i,j,k;ts.triangle_mesh.elements(t).Get(i,j,k);
        subdivided_values(ts.start_index_for_new_nodes-1+t)=(T)one_third*(base_values(i)+base_values(j)+base_values(k));}
}
}
}
#endif
