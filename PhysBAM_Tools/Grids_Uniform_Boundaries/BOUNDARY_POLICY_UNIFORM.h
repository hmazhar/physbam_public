//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BOUNDARY_POLICY_UNIFORM__
#define __BOUNDARY_POLICY_UNIFORM__

namespace PhysBAM{

template<class T_GRID,class T2> class BOUNDARY_UNIFORM;
template<class T_GRID,class T2> class BOUNDARY_REFLECTION_UNIFORM;
template<class TV> class GRID;
template<class T,int d> class VECTOR;

template<class T_GRID> struct BOUNDARY_POLICY;

template<class TV>
struct BOUNDARY_POLICY<GRID<TV> >
{
    typedef BOUNDARY_UNIFORM<GRID<TV>,typename TV::SCALAR> BOUNDARY_SCALAR;
    typedef BOUNDARY_REFLECTION_UNIFORM<GRID<TV>,typename TV::SCALAR> BOUNDARY_REFLECTION;
};

template<class T_BOUNDARY,class T_NEW> struct REBIND;
template<class T_GRID,class T2,class T_NEW> struct REBIND<BOUNDARY_UNIFORM<T_GRID,T2>,T_NEW>{typedef BOUNDARY_UNIFORM<T_GRID,T_NEW> TYPE;};
template<class T_BOUNDARY,int d_new> struct REBIND_LENGTH;
//template<class T_GRID,class T2,int d,int d_new> struct REBIND_LENGTH<BOUNDARY_UNIFORM<T_GRID,T2,d>,d_new>{typedef BOUNDARY_UNIFORM<T_GRID,T2,d_new> TYPE;};
template<class T_GRID,class T2,int d_new> struct REBIND_LENGTH<BOUNDARY_UNIFORM<T_GRID,T2>,d_new>{typedef BOUNDARY_UNIFORM<T_GRID,VECTOR<T2,d_new> > TYPE;};

//#####################################################################
}
#endif
