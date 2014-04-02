//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM
//#####################################################################
#ifndef __THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM__
#define __THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> //  T_AVERAGING=AVERAGING_UNIFORM<T_GRID>, T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2>
class THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM:public ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:
    template<class T3> struct REBIND{typedef THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T3,T_AVERAGING,typename T_INTERPOLATION::template REBIND<T3>::TYPE> TYPE;};
    template<class T_INTERPOLATION_2> struct REBIND_INTERPOLATION{typedef THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION_2> TYPE;};

    THREAD_QUEUE* thread_queue;
    int row_jump;

//#####################################################################
    THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM();
#ifdef USE_PTHREADS
    void Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost=0,const T_ARRAYS_T2* Z_max_ghost=0,T_ARRAYS_T2* Z_min=0,T_ARRAYS_T2* Z_max=0);
    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max);
     void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max);
//#####################################################################
#endif
};
}
#endif
