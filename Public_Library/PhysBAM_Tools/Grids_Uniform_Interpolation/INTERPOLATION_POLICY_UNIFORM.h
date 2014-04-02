//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERPOLATION_POLICY_UNIFORM__
#define __INTERPOLATION_POLICY_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> struct INTERPOLATION_POLICY;

template<class TV> class GRID;
template<class T,int d> class VECTOR;

template<class T>
struct INTERPOLATION_POLICY<GRID<VECTOR<T,1> > >
{
    typedef INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T> INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T> LINEAR_INTERPOLATION_SCALAR;
    typedef FACE_LOOKUP_UNIFORM<GRID<VECTOR<T,1> > > FACE_LOOKUP;
    typedef LINEAR_INTERPOLATION_MAC_1D_HELPER<GRID<VECTOR<T,1> > > LINEAR_INTERPOLATION_MAC_HELPER;
};

template<class T>
struct INTERPOLATION_POLICY<GRID<VECTOR<T,2> > >
{
    typedef INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T> INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T> LINEAR_INTERPOLATION_SCALAR;
    typedef FACE_LOOKUP_UNIFORM<GRID<VECTOR<T,2> > > FACE_LOOKUP;
    typedef LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<T,2> > > LINEAR_INTERPOLATION_MAC_HELPER;
};

template<class T>
struct INTERPOLATION_POLICY<GRID<VECTOR<T,3> > >
{
    typedef INTERPOLATION_UNIFORM<GRID<VECTOR<T,3> >,T> INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,3> >,T> LINEAR_INTERPOLATION_SCALAR;
    typedef FACE_LOOKUP_UNIFORM<GRID<VECTOR<T,3> > > FACE_LOOKUP;
    typedef LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<T,3> > > LINEAR_INTERPOLATION_MAC_HELPER;
};
}
#endif
