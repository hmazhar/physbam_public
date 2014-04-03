//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_DYADIC
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __VORTICITY_DYADIC__
#define __VORTICITY_DYADIC__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;

template<class T>
class VORTICITY_DYADIC
{
public:
//#####################################################################
    static void Vorticity(QUADTREE_GRID<T>& grid,const ARRAY<VECTOR<T,2> >& V_ghost,ARRAY<T>& vorticity);
    static void Vorticity(OCTREE_GRID<T>& grid,const ARRAY<VECTOR<T,3> >& V_ghost,ARRAY<VECTOR<T,3> >& vorticity);
//#####################################################################
};
}
#endif
#endif
