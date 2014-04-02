//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace HEAVISIDE
//##################################################################### 
#ifndef TRIANGULATED_SURFACE_HEAVISIDE_UNIFORM
#define TRIANGULATED_SURFACE_HEAVISIDE_UNIFORM
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{
template<class TV> class GRID;

namespace HEAVISIDE{
template<class T> void Calculate(TRIANGULATED_SURFACE<T>& surface,const GRID<VECTOR<T,3> >& grid,ARRAY<T,VECTOR<int,3> >& phi,bool print_progress=false);
};
};
#endif
