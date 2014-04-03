//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __ADVECTION_POLICY_RLE__
#define __ADVECTION_POLICY_RLE__

#include <PhysBAM_Tools/Advection/ADVECTION_FORWARD.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_RLE_FORWARD.h>

namespace PhysBAM{

template<class T_GRID>
struct ADVECTION_POLICY;

template<class T>
struct ADVECTION_POLICY<RLE_GRID_2D<T> >
{
private:
    typedef RLE_GRID_2D<T> T_GRID;
public:
    // normal
    typedef ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T> ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T> ADVECTION_HAMILTON_JACOBI_WENO_SCALAR;  // TODO: make this not completely wrong
};

template<class T>
struct ADVECTION_POLICY<RLE_GRID_3D<T> >
{
private:
    typedef RLE_GRID_3D<T> T_GRID;
public:
    // normal
    typedef ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T> ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T> ADVECTION_HAMILTON_JACOBI_WENO_SCALAR;  // TODO: make this not completely wrong
};

}
#endif
#endif
