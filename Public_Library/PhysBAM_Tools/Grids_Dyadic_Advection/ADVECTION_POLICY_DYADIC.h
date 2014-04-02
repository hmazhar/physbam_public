//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __ADVECTION_POLICY_DYADIC__
#define __ADVECTION_POLICY_DYADIC__

#include <PhysBAM_Tools/Advection/ADVECTION_FORWARD.h>
#include <PhysBAM_Tools/Grids_Dyadic_Advection/ADVECTION_DYADIC_FORWARD.h>

namespace PhysBAM{

template<class T_GRID>
struct ADVECTION_POLICY;

template<class T>
struct ADVECTION_POLICY<OCTREE_GRID<T> >
{
private:
    typedef OCTREE_GRID<T> T_GRID;
public:
    // normal
    typedef ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T> ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T> ADVECTION_HAMILTON_JACOBI_WENO_SCALAR; // TODO: make this not completely wrong
};

template<class T>
struct ADVECTION_POLICY<QUADTREE_GRID<T> >
{
private:
    typedef QUADTREE_GRID<T> T_GRID;
public:
    typedef ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T> ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T> ADVECTION_HAMILTON_JACOBI_WENO_SCALAR; // TODO: make this not completely wrong
};

}
#endif
#endif
