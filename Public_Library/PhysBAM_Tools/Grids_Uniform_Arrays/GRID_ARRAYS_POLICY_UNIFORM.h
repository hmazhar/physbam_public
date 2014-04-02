//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GRID_ARRAYS_POLICY_UNIFORM__
#define __GRID_ARRAYS_POLICY_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID> struct GRID_ARRAYS_POLICY;
template<class T,class ID> struct FACE_ARRAY;

template<class T>
struct GRID_ARRAYS_POLICY<GRID<VECTOR<T,1> > >
{
    typedef FLOOD_FILL_1D FLOOD_FILL;
    typedef ARRAY<T,VECTOR<int,1> > ARRAYS_SCALAR;
    typedef PhysBAM::ARRAY<T,FACE_INDEX<1> > FACE_ARRAYS;
    typedef ARRAYS_ND_BASE<VECTOR<T,1> > ARRAYS_BASE;
//#####################################################################
};

template<class T>
struct GRID_ARRAYS_POLICY<GRID<VECTOR<T,2> > >
{
    typedef FLOOD_FILL_2D FLOOD_FILL;
    typedef ARRAY<T,VECTOR<int,2> > ARRAYS_SCALAR;
    typedef PhysBAM::ARRAY<T,FACE_INDEX<2> > FACE_ARRAYS;
    typedef ARRAYS_ND_BASE<VECTOR<T,2> > ARRAYS_BASE;
//#####################################################################
};

template<class T>
struct GRID_ARRAYS_POLICY<GRID<VECTOR<T,3> > >
{
    typedef FLOOD_FILL_3D FLOOD_FILL;
    typedef ARRAY<T,VECTOR<int,3> > ARRAYS_SCALAR;
    typedef PhysBAM::ARRAY<T,FACE_INDEX<3> > FACE_ARRAYS;
    typedef ARRAYS_ND_BASE<VECTOR<T,3> > ARRAYS_BASE;
//#####################################################################
};

}
#endif
