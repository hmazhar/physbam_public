//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SEGMENTED_CURVE_INSIDE__
#define __SEGMENTED_CURVE_INSIDE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T,int d>
bool Inside_Any_Simplex(const SEGMENTED_CURVE<VECTOR<T,d> >& sc,const VECTOR<T,d>& location,int& segment_id,const T thickness_over_two)
{
    typedef typename BASIC_GEOMETRY_POLICY<VECTOR<T,d> >::SEGMENT T_SEGMENT;
    assert(sc.hierarchy);assert(sc.segment_list);
    ARRAY<int> nearby_segments;sc.hierarchy->Intersection_List(location,nearby_segments,thickness_over_two);
    for(int k=1;k<=nearby_segments.m;k++){
        T_SEGMENT& segment=(*sc.segment_list)(nearby_segments(k));
        if(segment.Inside(location,thickness_over_two)){segment_id=nearby_segments(k);return true;}}
    return false;
}
}
}
#endif
