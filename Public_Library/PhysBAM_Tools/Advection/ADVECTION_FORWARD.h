//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_FORWARD
//#####################################################################
#ifndef __ADVECTION_FORWARD__
#define __ADVECTION_FORWARD__

#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP=typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP> class ADVECTION;

}
#endif
