//#####################################################################
// Copyright 2009, Geoffrey Irving, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD
//#####################################################################
#ifndef __INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD__
#define __INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_FORWARD.h>

namespace PhysBAM{

template<class T_GRID,class T_NESTED_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID> > class FACE_LOOKUP_COLLIDABLE_UNIFORM;
template<class T_GRID,class T_NESTED_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID> > class FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM;
template<class T_GRID,class T_NESTED_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID> > class FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM;
template<class T_GRID,class T_FACE_LOOKUP> class AVERAGING_COLLIDABLE_UNIFORM;
template<class T_GRID,class T2> class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM;
template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> > class LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM;

}
#endif
