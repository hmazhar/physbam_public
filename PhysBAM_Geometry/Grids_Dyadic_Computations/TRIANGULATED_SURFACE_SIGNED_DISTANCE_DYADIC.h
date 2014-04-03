//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __TRIANGULATED_SURFACE_SIGNED_DISTANCE_DYADIC__
#define __TRIANGULATED_SURFACE_SIGNED_DISTANCE_DYADIC__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class T> class OCTREE_GRID;

namespace SIGNED_DISTANCE{
template<class T> void Calculate(TRIANGULATED_SURFACE<T>& surface,const GRID<VECTOR<T,3> >& uniform_grid,OCTREE_GRID<T>& grid,ARRAY<T>& phi,const int maximum_depth,const T half_bandwidth,bool verbose=false);
};
};
#endif
#endif
