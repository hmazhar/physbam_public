//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_STEP_PARTICLES
//#####################################################################
#ifndef __EULER_STEP_PARTICLES__
#define __EULER_STEP_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID>
class EULER_STEP_PARTICLES
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;
    typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR,TV>::TYPE T_ARRAYS_TV;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;
public:
//#####################################################################
    static void Euler_Step_Node(ARRAY_VIEW<TV> X,const T_GRID& grid,const T_ARRAYS_TV& U,const T dt);
    static void Euler_Step_Face(ARRAY_VIEW<TV> X,const T_GRID& grid,const T_FACE_ARRAYS& face_velocities,const T dt);
    static void Second_Order_Runge_Kutta_Step(ARRAY_VIEW<TV> X,const T_GRID& grid,const T_ARRAYS_TV& U,const T dt);
//#####################################################################
};
}
#endif
