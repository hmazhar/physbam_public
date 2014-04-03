//#####################################################################
// Copyright 2009, Geoffrey Irving, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header INTERPOLATION_COLLIDABLE_DYADIC_FORWARD
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __INTERPOLATION_COLLIDABLE_DYADIC_FORWARD__
#define __INTERPOLATION_COLLIDABLE_DYADIC_FORWARD__

#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC_FORWARD.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_FORWARD.h>

namespace PhysBAM{

template<class T_GRID,class T_NESTED_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID> > class FACE_LOOKUP_COLLIDABLE_DYADIC;
template<class T_GRID,class T_FACE_LOOKUP> class AVERAGING_COLLIDABLE_DYADIC;
template<class T_GRID,class T2> class LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC;
template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID> > class LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC;
}
#endif
#endif
