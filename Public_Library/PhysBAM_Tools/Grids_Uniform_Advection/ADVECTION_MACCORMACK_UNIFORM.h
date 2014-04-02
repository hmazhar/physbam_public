//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_MACCORMACK_UNIFORM
//#####################################################################
#ifndef __ADVECTION_MACCORMACK_UNIFORM__
#define __ADVECTION_MACCORMACK_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_NESTED_ADVECTION>
class ADVECTION_MACCORMACK_UNIFORM:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
private:
    // false to use 1st order advection
    const T_ARRAYS_BOOL* node_mask;
    const T_ARRAYS_BOOL* cell_mask;
    const T_FACE_ARRAYS_BOOL* face_mask;
    T_NESTED_ADVECTION& nested_advection;
public:
    bool clamp_extrema; // otherwise fall back to 1st order
    bool ensure_second_order; // sacrifice stability for accuracy.
    THREAD_QUEUE* thread_queue;

    ADVECTION_MACCORMACK_UNIFORM(T_NESTED_ADVECTION& nested_advection_input,const T_ARRAYS_BOOL* node_mask_input,const T_ARRAYS_BOOL* cell_mask_input,const T_FACE_ARRAYS_BOOL* face_mask_input,THREAD_QUEUE* thread_queue=0);

    void Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost_input=0,const T_ARRAYS_T2* Z_max_ghost_input=0,T_ARRAYS_T2* Z_min_input=0,T_ARRAYS_T2* Z_max_input=0);
    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const FACE_LOOKUP_UNIFORM<T_GRID>& face_velocities,
        T_BOUNDARY_T2& boundary,const T dt,const T time,const T_ARRAYS_T2* Z_min_ghost_input,const T_ARRAYS_T2* Z_max_ghost_input,T_ARRAYS_T2* Z_min_input,
        T_ARRAYS_T2* Z_max_input);
    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const FACE_LOOKUP_UNIFORM<T_GRID>& Z_ghost,
        const FACE_LOOKUP_UNIFORM<T_GRID>& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,const FACE_LOOKUP_UNIFORM<T_GRID>* Z_min_ghost_input,
        const FACE_LOOKUP_UNIFORM<T_GRID>* Z_max_ghost_input,T_FACE_ARRAYS_SCALAR* Z_min_input,T_FACE_ARRAYS_SCALAR* Z_max_input);
private:
    void Apply_Clamped_Extrema_Limiter_Node_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_forward_ghost,const T_ARRAYS_T2& Z_backward_ghost,const T_ARRAYS_T2& Z_min,const T_ARRAYS_T2& Z_max);
    void Apply_Clamped_Extrema_Limiter_Cell_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_forward_ghost,const T_ARRAYS_T2& Z_backward_ghost,const T_ARRAYS_T2& Z_min,const T_ARRAYS_T2& Z_max);
    void Apply_Clamped_Extrema_Limiter_Face_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost,const T_FACE_ARRAYS_SCALAR& Z_min,const T_FACE_ARRAYS_SCALAR& Z_max);

    void Apply_Reversion_Limiter_Node_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_forward_ghost,const T_ARRAYS_T2& Z_backward_ghost,const T_ARRAYS_T2& Z_min,const T_ARRAYS_T2& Z_max);
    void Apply_Reversion_Limiter_Cell_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_forward_ghost,const T_ARRAYS_T2& Z_backward_ghost,const T_ARRAYS_T2& Z_min,const T_ARRAYS_T2& Z_max);
    void Apply_Reversion_Limiter_Face_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost,const T_FACE_ARRAYS_SCALAR& Z_min,const T_FACE_ARRAYS_SCALAR& Z_max);

    void Apply_Second_Order_Update_Face_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost);

//#####################################################################
};
}
#endif
