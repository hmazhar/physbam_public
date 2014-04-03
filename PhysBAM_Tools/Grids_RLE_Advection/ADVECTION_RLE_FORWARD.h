//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_RLE_FORWARD
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __ADVECTION_RLE_FORWARD__
#define __ADVECTION_RLE_FORWARD__

#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_RLE<T_GRID>,class T_INTERPOLATION=LINEAR_INTERPOLATION_RLE<T_GRID,T2> > class ADVECTION_SEMI_LAGRANGIAN_RLE;

template<class T> class ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE;

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION,class T_NEW> struct REBIND<ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>,T_NEW>{
    typedef ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T_NEW,T_AVERAGING,typename REBIND<T_INTERPOLATION,T_NEW>::TYPE> TYPE;};

}
#endif
#endif
