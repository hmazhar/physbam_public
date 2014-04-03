//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_POLICY 
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_POLICY__
#define __RLE_GRID_POLICY__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
namespace PhysBAM{

template<class TV> struct RLE_TAG{};

template<class T> class RLE_GRID_1D;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

//#####################################################################

template<class TV> struct RLE_GRID_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct RLE_GRID_POLICY<VECTOR<T,0> >
{
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct RLE_GRID_POLICY<VECTOR<T,1> >
{
    typedef RLE_GRID_1D<T> RLE_GRID;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct RLE_GRID_POLICY<VECTOR<T,2> >
{
    typedef RLE_GRID_2D<T> RLE_GRID;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct RLE_GRID_POLICY<VECTOR<T,3> >
{
    typedef RLE_GRID_3D<T> RLE_GRID;
};
//#####################################################################

}
#endif
#endif
