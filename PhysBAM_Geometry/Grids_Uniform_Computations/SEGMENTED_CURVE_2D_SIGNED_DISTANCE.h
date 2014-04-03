//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef SEGMENTED_CURVE_2D_SIGNED_DISTANCE
#define SEGMENTED_CURVE_2D_SIGNED_DISTANCE
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{
template<class TV> class GRID;

namespace SIGNED_DISTANCE{
template<class T> void Calculate(SEGMENTED_CURVE_2D<T>& curve,const GRID<VECTOR<T,2> >& grid,ARRAY<T,VECTOR<int,2> >& phi,bool print_progress=false);
};
};
#endif
