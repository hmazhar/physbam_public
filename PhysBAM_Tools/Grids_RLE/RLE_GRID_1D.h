//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_1D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_1D__
#define __RLE_GRID_1D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID.h>

namespace PhysBAM{

template<class T>
class RLE_GRID_1D
{
public:
    typedef int INDEX;
    typedef T SCALAR;
    typedef VECTOR<T,1> VECTOR_T;
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef RLE_TAG<VECTOR_T> GRID_TAG;

    RLE_GRID_1D();
    ~RLE_GRID_1D();
};
}
#endif
#endif
