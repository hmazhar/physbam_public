//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_COLLIDABLE_RLE_FORWARD
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __ADVECTION_COLLIDABLE_RLE_FORWARD__
#define __ADVECTION_COLLIDABLE_RLE_FORWARD__

#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_RLE_FORWARD.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_RLE<T_GRID> > class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE;
template<class T_GRID,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_RLE<T_GRID> > class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE;
}
#endif
#endif
