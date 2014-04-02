//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADVECTION_POLICY_UNIFORM__
#define __ADVECTION_POLICY_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class T_GRID>
struct ADVECTION_POLICY
{
private:
    typedef typename T_GRID::SCALAR T;
public:
    // normal
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<T_GRID,T> ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ::PhysBAM::ADVECTION_CONSERVATIVE_ENO<T_GRID,T> ADVECTION_CONSERVATIVE_ENO;
    typedef ADVECTION_HAMILTON_JACOBI_WENO<T_GRID,T> ADVECTION_HAMILTON_JACOBI_WENO_SCALAR;
};

}
#endif
