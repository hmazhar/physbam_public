//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING_METHOD_UNIFORM  
//#####################################################################
#ifndef __FAST_MARCHING_METHOD_UNIFORM__
#define __FAST_MARCHING_METHOD_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class FAST_MARCHING_METHOD_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
public:
    const T_LEVELSET& levelset;
private:
    T_GRID cell_grid;
    TV_INT dimension_start,dimension_end;
    int ghost_cells;
public:
    THREAD_QUEUE* thread_queue;

    FAST_MARCHING_METHOD_UNIFORM(const T_LEVELSET& levelset,const int ghost_cells,THREAD_QUEUE* thread_queue_input=0);
    ~FAST_MARCHING_METHOD_UNIFORM();

    bool Neighbor_Visible(const int neighbor_number,const TV_INT& current_index) // neighbor_number between 1 and 3 -- right, top, back
    {return !levelset.collision_body_list->cell_neighbors_visible.Valid_Index(current_index) || levelset.collision_body_list->cell_neighbors_visible(current_index)(neighbor_number);}

//#####################################################################
    void Fast_Marching_Method_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells=false);
    void Fast_Marching_Method(T_ARRAYS_SCALAR& phi_ghost,const T stopping_distance=0,const ARRAY<TV_INT>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Fast_Marching_Method(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& seed_indices,const T stopping_distance=0,const bool add_seed_indices_for_ghost_cells=false);
    void Initialize_Interface_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_SCALAR& phi_new,T_ARRAYS_BOOL& done);
private:
    void Update_Or_Add_Neighbor(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const TV_INT& neighbor);
    void Initialize_Interface(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells=false);
    void Initialize_Interface(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices=0,
        const bool add_seed_indices_for_ghost_cells=false);
    void Initialize_Interface(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const bool add_seed_indices_for_ghost_cells=false);
    void Update_Close_Point(T_ARRAYS_SCALAR& phi_ghost,const T_ARRAYS_BOOL& done,const TV_INT& index);
    void Add_To_Initial(T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,const TV_INT& index);
//#####################################################################
};
}
#endif
