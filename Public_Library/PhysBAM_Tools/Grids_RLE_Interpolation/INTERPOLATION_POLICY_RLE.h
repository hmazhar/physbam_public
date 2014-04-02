//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __INTERPOLATION_POLICY_RLE__
#define __INTERPOLATION_POLICY_RLE__

#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_RLE_FORWARD.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> struct INTERPOLATION_POLICY;

template<class TV> class GRID;
template<class T> class RLE_GRID_1D;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;
template<class T,int d> class VECTOR;

template<class T>
struct INTERPOLATION_POLICY<RLE_GRID_2D<T> >
{
private:
    typedef RLE_GRID_2D<T> T_GRID;
public:
    typedef INTERPOLATION_RLE<T_GRID,T> INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_RLE_HELPER<T_GRID> LINEAR_INTERPOLATION_HELPER;
    typedef LINEAR_INTERPOLATION_MAC_2D_HELPER<T_GRID> LINEAR_INTERPOLATION_MAC_HELPER;
    typedef FACE_LOOKUP_RLE<T_GRID> FACE_LOOKUP;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T> LINEAR_INTERPOLATION_SCALAR;
};

template<class T>
struct INTERPOLATION_POLICY<RLE_GRID_3D<T> >
{
private:
    typedef RLE_GRID_3D<T> T_GRID;
public:
    typedef INTERPOLATION_RLE<T_GRID,T> INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_RLE_HELPER<T_GRID> LINEAR_INTERPOLATION_HELPER;
    typedef LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID> LINEAR_INTERPOLATION_MAC_HELPER;
    typedef FACE_LOOKUP_RLE<T_GRID> FACE_LOOKUP;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T> LINEAR_INTERPOLATION_SCALAR;
};
}
#endif
#endif
