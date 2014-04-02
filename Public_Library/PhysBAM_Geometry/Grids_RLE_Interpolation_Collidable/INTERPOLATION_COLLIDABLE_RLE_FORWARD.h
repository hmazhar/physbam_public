//#####################################################################
// Copyright 2009, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header INTERPOLATION_COLLIDABLE_RLE_FORWARD
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __INTERPOLATION_COLLIDABLE_RLE_FORWARD__
#define __INTERPOLATION_COLLIDABLE_RLE_FORWARD__

#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_RLE_FORWARD.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_FORWARD.h>

namespace PhysBAM{

template<class T_GRID,class T_NESTED_LOOKUP=FACE_LOOKUP_RLE<T_GRID> > class FACE_LOOKUP_COLLIDABLE_RLE;
template<class T_GRID,class T_FACE_LOOKUP> class AVERAGING_COLLIDABLE_RLE;
template<class T_GRID,class T2> class LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE;
template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_RLE<T_GRID> > class LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE;

}
#endif
#endif
