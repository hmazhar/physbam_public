//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_ON_A_RAY_POLICY 
//#####################################################################
#ifndef __IMPLICIT_OBJECT_ON_A_RAY_POLICY__
#define __IMPLICIT_OBJECT_ON_A_RAY_POLICY__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_FORWARD.h>
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_POLICY.h>
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_POLICY.h>
#endif
namespace PhysBAM{

template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T_IMPLICIT_OBJECT> class IMPLICIT_OBJECT_ON_A_RAY;

template<class T_GRID> struct IMPLICIT_OBJECT_ON_A_RAY_POLICY:public IMPLICIT_OBJECT_ON_A_RAY_POLICY<typename T_GRID::GRID_TAG>{};

//#####################################################################
// Uniform
//#####################################################################
template<class TV>
struct IMPLICIT_OBJECT_ON_A_RAY_POLICY<UNIFORM_TAG<TV> >
{
    typedef IMPLICIT_OBJECT_ON_A_RAY<LEVELSET_IMPLICIT_OBJECT<TV> > TYPE;
};
//#####################################################################
// Dyadic
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class TV> class DYADIC_IMPLICIT_OBJECT_ON_A_RAY;
template<class TV>
struct IMPLICIT_OBJECT_ON_A_RAY_POLICY<DYADIC_TAG<TV> >
{
    typedef DYADIC_IMPLICIT_OBJECT_ON_A_RAY<typename DYADIC_GRID_POLICY<TV>::DYADIC_GRID> TYPE;
};
#endif
//#####################################################################
// RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class TV> class RLE_LEVELSET_ON_A_RAY;
template<class TV>
struct IMPLICIT_OBJECT_ON_A_RAY_POLICY<RLE_TAG<TV> >
{
    typedef RLE_LEVELSET_ON_A_RAY<TV> TYPE;
};
#endif
//#####################################################################
}
#endif
