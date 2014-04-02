//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_ADVECTION_SEMI_LAGRANGIAN_TASK.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

using namespace PhysBAM;

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM()
    :thread_queue(0),row_jump(1)
{}

#ifdef USE_PTHREADS
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    RANGE<TV_INT> domain(grid.Domain_Indices());
    int min_value=domain.min_corner.x,max_value=domain.max_corner.x;
    for(int i=min_value;i<=max_value;i+=row_jump){
        domain.min_corner.x=min(i,max_value);domain.max_corner.x=min(i+row_jump-1,max_value);
        ADVECTION_SEMI_LAGRANGIAN_TASK_NODE<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>* task=
            new ADVECTION_SEMI_LAGRANGIAN_TASK_NODE<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>(grid,Z,Z_ghost,V,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max,domain);
        thread_queue->Queue(task);}
    thread_queue->Wait();
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    RANGE<TV_INT> domain(grid.Domain_Indices());
    int min_value=domain.min_corner.x,max_value=domain.max_corner.x;
    for(int i=min_value;i<=max_value;i+=row_jump){
        domain.min_corner.x=min(i,max_value);domain.max_corner.x=min(i+row_jump-1,max_value);
        ADVECTION_SEMI_LAGRANGIAN_TASK_CELL<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>* task=
            new ADVECTION_SEMI_LAGRANGIAN_TASK_CELL<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>(grid,Z,Z_ghost,face_velocities,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max,domain);
        thread_queue->Queue(task);}
    thread_queue->Wait();
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    for(int i=1;i<=TV::dimension;i++){
        RANGE<TV_INT> domain(grid.Domain_Indices());
        int min_value=domain.min_corner(i),max_value=domain.max_corner(i)+1;
        for(int j=min_value;j<=max_value;j+=row_jump){
            domain.min_corner(i)=j;domain.max_corner(i)=min(j+row_jump-1,max_value);
            ADVECTION_SEMI_LAGRANGIAN_TASK_FACE<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>* task=
                new ADVECTION_SEMI_LAGRANGIAN_TASK_FACE<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>(grid,Z,Z_ghost,face_velocities,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max,domain,i);
            thread_queue->Queue(task);}}
    thread_queue->Wait();
}
#endif
