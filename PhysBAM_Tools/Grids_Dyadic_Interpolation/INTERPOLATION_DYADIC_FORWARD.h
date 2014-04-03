//#####################################################################
// Copyright 2009, Andrew Selle, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __INTERPOLATION_DYADIC_FORWARD__
#define __INTERPOLATION_DYADIC_FORWARD__

namespace PhysBAM{
template<class T_GRID> class FACE_LOOKUP_DYADIC;

template<class T_GRID> class LINEAR_INTERPOLATION_MAC_1D_HELPER;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_2D_HELPER;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_3D_HELPER;
template<class T_GRID> class LINEAR_INTERPOLATION_DYADIC_HELPER;
template<class T> class LINEAR_INTERPOLATION_BINTREE_HELPER;
template<class T> class LINEAR_INTERPOLATION_OCTREE_HELPER;
template<class T> class LINEAR_INTERPOLATION_QUADTREE_HELPER;


template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID> > class INTERPOLATION_DYADIC;
template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID> > class LINEAR_INTERPOLATION_DYADIC;

template<class T_GRID,class T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID> > class AVERAGING_DYADIC;

template<class T_INTERPOLATION,class T_NEW> struct REBIND;
template<class T_GRID,class T2,class T_LOOKUP,class T_NEW> struct REBIND<INTERPOLATION_DYADIC<T_GRID,T2,T_LOOKUP>,T_NEW>{typedef INTERPOLATION_DYADIC<T_GRID,T_NEW,T_LOOKUP> TYPE;};
template<class T_GRID,class T2,class T_LOOKUP,class T_NEW> struct REBIND<LINEAR_INTERPOLATION_DYADIC<T_GRID,T2,T_LOOKUP>,T_NEW>{typedef LINEAR_INTERPOLATION_DYADIC<T_GRID,T_NEW,T_LOOKUP> TYPE;};

template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID> > class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_DYADIC;

}
#endif
#endif
