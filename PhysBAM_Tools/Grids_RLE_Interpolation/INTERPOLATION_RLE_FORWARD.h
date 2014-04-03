//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __INTERPOLATION_RLE_FORWARD__
#define __INTERPOLATION_RLE_FORWARD__

namespace PhysBAM{
template<class T_GRID> class FACE_LOOKUP_RLE;

template<class T_GRID> class LINEAR_INTERPOLATION_MAC_1D_HELPER;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_2D_HELPER;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_3D_HELPER;
template<class T_GRID> class LINEAR_INTERPOLATION_RLE_HELPER;

template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_RLE<T_GRID> > class INTERPOLATION_RLE;
template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_RLE<T_GRID> > class LINEAR_INTERPOLATION_RLE;

template<class T_INTERPOLATION,class T_NEW> struct REBIND;
template<class T_GRID,class T2,class T_LOOKUP,class T_NEW> struct REBIND<INTERPOLATION_RLE<T_GRID,T2,T_LOOKUP>,T_NEW>{typedef INTERPOLATION_RLE<T_GRID,T_NEW,T_LOOKUP> TYPE;};
template<class T_GRID,class T2,class T_LOOKUP,class T_NEW> struct REBIND<LINEAR_INTERPOLATION_RLE<T_GRID,T2,T_LOOKUP>,T_NEW>{typedef LINEAR_INTERPOLATION_RLE<T_GRID,T_NEW,T_LOOKUP> TYPE;};

template<class T_GRID,class T_FACE_LOOKUP=FACE_LOOKUP_RLE<T_GRID> > class AVERAGING_RLE;

}
#endif
#endif
