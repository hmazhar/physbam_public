//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __BOUNDARY_POLICY_DYADIC__
#define __BOUNDARY_POLICY_DYADIC__

namespace PhysBAM{

template<class T_GRID,class T2> class BOUNDARY_DYADIC;
template<class T> class BINTREE_GRID;

template<class T_GRID> struct BOUNDARY_POLICY;

template<class T>
struct BOUNDARY_POLICY<BINTREE_GRID<T> >
{
    typedef BOUNDARY_DYADIC<BINTREE_GRID<T>,T> BOUNDARY_SCALAR;
    typedef BOUNDARY_DYADIC<BINTREE_GRID<T>,T> BOUNDARY_REFLECTION; // TODO: actually set this to a reflecting boundary
};


#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;

template<class T>
struct BOUNDARY_POLICY<QUADTREE_GRID<T> >
{
    typedef BOUNDARY_DYADIC<QUADTREE_GRID<T>,T> BOUNDARY_SCALAR;
    typedef BOUNDARY_DYADIC<QUADTREE_GRID<T>,T> BOUNDARY_REFLECTION; // TODO: actually set this to a reflecting boundary
};

template<class T_GRID> struct BOUNDARY_POLICY;

template<class T>
struct BOUNDARY_POLICY<OCTREE_GRID<T> >
{
    typedef BOUNDARY_DYADIC<OCTREE_GRID<T>,T> BOUNDARY_SCALAR;
    typedef BOUNDARY_DYADIC<OCTREE_GRID<T>,T> BOUNDARY_REFLECTION; // TODO: actually set this to a reflecting boundary
};
#endif

template<class T_BOUNDARY,class T_NEW> struct REBIND;
template<class T_GRID,class T2,class T_NEW> struct REBIND<BOUNDARY_DYADIC<T_GRID,T2>,T_NEW>{typedef BOUNDARY_DYADIC<T_GRID,T_NEW> TYPE;};
template<class T_BOUNDARY,int d_new> struct REBIND_LENGTH;
template<class T_GRID,class T2,int d_new> struct REBIND_LENGTH<BOUNDARY_DYADIC<T_GRID,T2>,d_new>{typedef BOUNDARY_DYADIC<T_GRID,T2> TYPE;};

//#####################################################################
}
#endif
#endif
