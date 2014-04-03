//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TETRAHEDRALIZED_VOLUME_INSIDE__
#define __TETRAHEDRALIZED_VOLUME_INSIDE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T>
bool Inside_Any_Simplex(const TETRAHEDRALIZED_VOLUME<T>& tv,const VECTOR<T,3>& location,int& tetrahedron_id,const T thickness_over_two)
{
    assert(tv.hierarchy);assert(tv.tetrahedron_list);
    ARRAY<int> nearby_tetrahedrons;tv.hierarchy->Intersection_List(location,nearby_tetrahedrons,thickness_over_two);
    for(int k=1;k<=nearby_tetrahedrons.m;k++){
        TETRAHEDRON<T>& tetrahedron=(*tv.tetrahedron_list)(nearby_tetrahedrons(k));
        if(tetrahedron.Inside(location,thickness_over_two)){tetrahedron_id=nearby_tetrahedrons(k);return true;}}
    return false;
}
}
}
#endif
