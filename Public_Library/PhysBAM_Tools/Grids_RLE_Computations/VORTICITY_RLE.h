//#####################################################################
// Copyright 2009, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_RLE
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __VORTICITY_RLE__
#define __VORTICITY_RLE__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;

template<class TV>
class VORTICITY_RLE
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename RLE_GRID_POLICY<TV>::RLE_GRID RLE_GRID;
public:
//#####################################################################
    static void Vorticity(const RLE_GRID& grid,const ARRAY<T>& V,ARRAY<T_SPIN>& vorticity);
//#####################################################################
};
}
#endif
#endif
