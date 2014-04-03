//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_UNIFORM_FORWARD
//#####################################################################
#ifndef __ADVECTION_UNIFORM_FORWARD__
#define __ADVECTION_UNIFORM_FORWARD__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_UNIFORM<T_GRID>,class T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2> > class ADVECTION_SEMI_LAGRANGIAN_UNIFORM;
template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_UNIFORM<T_GRID>,class T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2> > class ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA;
template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_UNIFORM<T_GRID>,class T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2> > class THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM;

template<class T_GRID,class T2,class T_NESTED_ADVECTION> class ADVECTION_MACCORMACK_UNIFORM;

template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_UNIFORM<T_GRID> > class ADVECTION_SEPARABLE_UNIFORM;
template<class T,class T2> class ADVECTION_CONSERVATIVE_ENO;
template<class T,class T2> class ADVECTION_CONSERVATIVE_WENO;
template<class T_GRID,class T2> class ADVECTION_CENTRAL;
template<class T_GRID,class T2> class ADVECTION_HAMILTON_JACOBI_ENO;
template<class T_GRID,class T2> class ADVECTION_HAMILTON_JACOBI_WENO;
template<class T,class T2> class ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO;
template<class T,class T2> class ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO;
template<class T_GRID> class ADVECTION_IMPROVED_VOF_WITH_REFINEMENT;

template<class T_ADVECTION,class T_NEW> struct REBIND;
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION,class T_NEW> struct REBIND<ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>,T_NEW>{
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T_NEW,T_AVERAGING,typename REBIND<T_INTERPOLATION,T_NEW>::TYPE> TYPE;};
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION,class T_NEW> struct REBIND<ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>,T_NEW>{
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<T_GRID,T_NEW,T_AVERAGING,typename REBIND<T_INTERPOLATION,T_NEW>::TYPE> TYPE;};

}
#endif
