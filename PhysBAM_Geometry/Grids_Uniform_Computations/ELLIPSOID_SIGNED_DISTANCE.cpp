//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/ELLIPSOID_SIGNED_DISTANCE.h>

namespace PhysBAM{
namespace SIGNED_DISTANCE{
//#####################################################################
// Function Calculate
//#####################################################################
template<class TV,class T_GRID,class T_ARRAY> void Calculate_Approximate(const ELLIPSOID<TV>& ellipsoid,const T_GRID& grid,T_ARRAY& phi,bool verbose)
{
    for(int i=1;i<=grid.counts.x;i++)for(int j=1;j<=grid.counts.y;j++)for(int ij=1;ij<=grid.counts.z;ij++)
        phi(i,j,ij)=ellipsoid.Approximate_Signed_Distance(grid.X(i,j,ij));
}
//#####################################################################

template void Calculate_Approximate(const ELLIPSOID<float>&,const GRID<VECTOR<float,3> >&,GRID_ARRAYS_POLICY<GRID<VECTOR<float,3> > >::ARRAYS_SCALAR&,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Calculate_Approximate(const ELLIPSOID<double>&,const GRID<VECTOR<double,3> >&,GRID_ARRAYS_POLICY<GRID<VECTOR<double,3> > >::ARRAYS_SCALAR&,bool);
#endif
}
}
