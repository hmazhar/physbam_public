//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __LEVELSET_DYADIC_SIGNED_DISTANCE__
#define __LEVELSET_DYADIC_SIGNED_DISTANCE__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{
template<class T> class OCTREE_GRID;
template<class T_GRID> class LEVELSET_DYADIC;
template<class T_GRID> class LEVELSET_3D;

namespace SIGNED_DISTANCE{
template<class T_GRID,class T_UNIFORM,class T_ARRAYS> void Calculate(const LEVELSET_DYADIC<T_GRID>& levelset,const T_UNIFORM& grid,T_ARRAYS& phi,bool verbose=false);
template<class T_GRID,class T> void Calculate(const LEVELSET_3D<T_GRID>& levelset,OCTREE_GRID<T>& grid,ARRAY<T>& phi,bool verbose=false);
template<class T_GRID,class T> void Calculate_Based_On_Curvature(const LEVELSET_3D<T_GRID>& levelset,const OCTREE_GRID<T>& grid,ARRAY<T,VECTOR<int,1> >& phi,const int levels,bool output_stastics=false);
};
};
#endif
#endif
