//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __GRID_ARRAYS_POLICY__
#define __GRID_ARRAYS_POLICY__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM{

template<class T> class RLE_GRID_1D;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T>
struct GRID_ARRAYS_POLICY<RLE_GRID_1D<T> >
{
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef ARRAY<T> FACE_ARRAYS;
    typedef ARRAY_BASE<T,ARRAY<T,int>,int> ARRAYS_BASE;
//#####################################################################
};

template<class T>
struct GRID_ARRAYS_POLICY<RLE_GRID_2D<T> >
{
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef ARRAY<T> FACE_ARRAYS;
    typedef ARRAY_BASE<T,ARRAY<T,int>,int> ARRAYS_BASE;
//#####################################################################
};

template<class T>
struct GRID_ARRAYS_POLICY<RLE_GRID_3D<T> >
{
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef ARRAY<T> FACE_ARRAYS;
    typedef ARRAY_BASE<T,ARRAY<T,int>,int> ARRAYS_BASE;
//#####################################################################
};

/*template<class T>
struct GRID_ARRAYS_POLICY<RLE_GRID<T> >
{
    typedef ARRAY<T> ARRAYS_SCALAR;
    typedef ARRAY<T> FACE_ARRAYS;
//#####################################################################
};*/

}
#endif
#endif
