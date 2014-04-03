//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_DYADIC_FORWARD
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __ADVECTION_DYADIC_FORWARD__
#define __ADVECTION_DYADIC_FORWARD__

#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC_FORWARD.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_DYADIC<T_GRID>,class T_INTERPOLATION=LINEAR_INTERPOLATION_DYADIC<T_GRID,T2> > class ADVECTION_SEMI_LAGRANGIAN_DYADIC;
template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_DYADIC<T_GRID> > class ADVECTION_HAMILTON_JACOBI_ENO_DYADIC;

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION,class T_NEW> struct REBIND<ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>,T_NEW>{
    typedef ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T_NEW,T_AVERAGING,typename REBIND<T_INTERPOLATION,T_NEW>::TYPE> TYPE;};

}
#endif
#endif
