//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __BOUNDARY_POLICY_RLE__
#define __BOUNDARY_POLICY_RLE__

namespace PhysBAM{

template<class T_GRID,class T2> class BOUNDARY_RLE;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class T_GRID> struct BOUNDARY_POLICY;

template<class T>
struct BOUNDARY_POLICY<RLE_GRID_2D<T> >
{
    typedef BOUNDARY_RLE<RLE_GRID_2D<T>,T> BOUNDARY_SCALAR;
    typedef BOUNDARY_RLE<RLE_GRID_2D<T>,T> BOUNDARY_REFLECTION; // TODO: actually set this to a reflecting boundary
};

template<class T>
struct BOUNDARY_POLICY<RLE_GRID_3D<T> >
{
    typedef BOUNDARY_RLE<RLE_GRID_3D<T>,T> BOUNDARY_SCALAR;
    typedef BOUNDARY_RLE<RLE_GRID_3D<T>,T> BOUNDARY_REFLECTION; // TODO: actually set this to a reflecting boundary
};

template<class T_BOUNDARY,class T_NEW> struct REBIND;
template<class T_GRID,class T2,class T_NEW> struct REBIND<BOUNDARY_RLE<T_GRID,T2>,T_NEW>{typedef BOUNDARY_RLE<T_GRID,T_NEW> TYPE;};
template<class T_BOUNDARY,int d_new> struct REBIND_LENGTH;
template<class T_GRID,class T2,int d_new> struct REBIND_LENGTH<BOUNDARY_RLE<T_GRID,T2>,d_new>{typedef BOUNDARY_RLE<T_GRID,T2> TYPE;};
//#####################################################################
}
#endif
#endif
