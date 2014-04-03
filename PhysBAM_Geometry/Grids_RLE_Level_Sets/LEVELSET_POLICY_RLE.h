//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_POLICY_RLE 
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __LEVELSET_POLICY_RLE__
#define __LEVELSET_POLICY_RLE__

#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
namespace PhysBAM{

template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class TV> class RLE_IMPLICIT_OBJECT;

template<class T_GRID> class LEVELSET_RLE;
template<class T_GRID> class FAST_LEVELSET;
template<class T_GRID> class LEVELSET_MULTIPLE_RLE;
template<class T_GRID> class PARTICLE_LEVELSET_RLE;
template<class T_GRID> class PARTICLE_LEVELSET_EVOLUTION_RLE;
template<class T_GRID,class T2> struct EXTRAPOLATION_RLE;

template<class T> struct LEVELSET_POLICY<RLE_GRID_2D<T> >
{
    typedef LEVELSET_RLE<RLE_GRID_2D<T> > LEVELSET;
    typedef LEVELSET FAST_LEVELSET_T;
    typedef LEVELSET_MULTIPLE_RLE<RLE_GRID_2D<T> > LEVELSET_MULTIPLE;
    typedef RLE_IMPLICIT_OBJECT<VECTOR<T,2> > LEVELSET_IMPLICIT_OBJECT;
    typedef PARTICLE_LEVELSET_RLE<RLE_GRID_2D<T> > PARTICLE_LEVELSET;
    typedef PARTICLE_LEVELSET_EVOLUTION_RLE<RLE_GRID_2D<T> > PARTICLE_LEVELSET_EVOLUTION;
};

template<class T> struct LEVELSET_POLICY<RLE_GRID_3D<T> >
{
    typedef LEVELSET_RLE<RLE_GRID_3D<T> > LEVELSET;
    typedef LEVELSET FAST_LEVELSET_T;
    typedef LEVELSET_MULTIPLE_RLE<RLE_GRID_3D<T> > LEVELSET_MULTIPLE;
    typedef RLE_IMPLICIT_OBJECT<VECTOR<T,3> > LEVELSET_IMPLICIT_OBJECT;
    typedef PARTICLE_LEVELSET_RLE<RLE_GRID_3D<T> > PARTICLE_LEVELSET;
    typedef PARTICLE_LEVELSET_EVOLUTION_RLE<RLE_GRID_3D<T> > PARTICLE_LEVELSET_EVOLUTION;
};

}
#endif
#endif
