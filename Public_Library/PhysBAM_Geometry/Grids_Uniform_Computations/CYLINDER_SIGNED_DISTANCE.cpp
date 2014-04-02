//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/CYLINDER_SIGNED_DISTANCE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/SIGNED_DISTANCE.h>

namespace PhysBAM{

template<class TV> class GRID;

namespace SIGNED_DISTANCE{
#define CYLINDER_SIGNED_DISTANCE_HELPER(T) \
    template void Calculate(CYLINDER<T>&,const GRID<VECTOR<T,3> >&,GRID_ARRAYS_POLICY<GRID<VECTOR<T,3> > >::ARRAYS_SCALAR&,bool);

CYLINDER_SIGNED_DISTANCE_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
CYLINDER_SIGNED_DISTANCE_HELPER(double);
#endif
}
}
