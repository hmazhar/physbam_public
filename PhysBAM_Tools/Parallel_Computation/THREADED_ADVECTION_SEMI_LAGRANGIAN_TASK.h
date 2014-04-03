//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_PTHREADS
#ifndef __THREADED_ADVECTION_SEMI_LAGRANGIAN_TASK__
#define __THREADED_ADVECTION_SEMI_LAGRANGIAN_TASK__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION>
class ADVECTION_SEMI_LAGRANGIAN_TASK_NODE:public THREAD_QUEUE::TASK
{    
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;

public:

    const T_GRID& grid;
    T_ARRAYS_T2& Z;
    const T_ARRAYS_T2& Z_ghost;
    const T_ARRAYS_VECTOR& V;
    T_BOUNDARY_T2& boundary;
    const T dt;
    const T time;
    const T_ARRAYS_T2* Z_min_ghost;
    const T_ARRAYS_T2* Z_max_ghost;
    T_ARRAYS_T2* Z_min;
    T_ARRAYS_T2* Z_max;
    RANGE<TV_INT> domain;

    ADVECTION_SEMI_LAGRANGIAN_TASK_NODE(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max,RANGE<TV_INT> domain)
        :grid(grid),Z(Z),Z_ghost(Z_ghost),V(V),boundary(boundary),dt(dt),time(time),Z_min_ghost(Z_min_ghost),Z_max_ghost(Z_max_ghost),Z_min(Z_min),Z_max(Z_max),domain(domain)
    {}
    
    void Run()
    {
        T_INTERPOLATION interpolation;
        if(Z_min && Z_max) for(NODE_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
            TV_INT node=iterator.Node_Index();TV X=iterator.Location()-dt*V(node);
            Z(node)=interpolation.Clamped_To_Array_Node(grid,Z_ghost,X);
            VECTOR<T2,2> extrema=interpolation.Extrema_Clamped_To_Array_Node(grid,*Z_min_ghost,*Z_max_ghost,X);
            (*Z_min)(node)=extrema.x;(*Z_max)(node)=extrema.y;}
        else for(NODE_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            Z(node)=interpolation.Clamped_To_Array(grid,Z_ghost,iterator.Location()-dt*V(node));}
    }
};

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION>
class ADVECTION_SEMI_LAGRANGIAN_TASK_CELL:public THREAD_QUEUE::TASK
{    
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;

public:

    const T_GRID& grid;
    T_ARRAYS_T2& Z;
    const T_ARRAYS_T2& Z_ghost;
    const T_FACE_LOOKUP& face_velocities;
    T_BOUNDARY_T2& boundary;
    const T dt;
    const T time;
    const T_ARRAYS_T2* Z_min_ghost;
    const T_ARRAYS_T2* Z_max_ghost;
    T_ARRAYS_T2* Z_min;
    T_ARRAYS_T2* Z_max;
    RANGE<TV_INT> domain;

    ADVECTION_SEMI_LAGRANGIAN_TASK_CELL(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max,RANGE<TV_INT> domain)
        :grid(grid),Z(Z),Z_ghost(Z_ghost),face_velocities(face_velocities),boundary(boundary),dt(dt),time(time),Z_min_ghost(Z_min_ghost),Z_max_ghost(Z_max_ghost),Z_min(Z_min),Z_max(Z_max),domain(domain)
    {}
    
    void Run()
    {
        T_INTERPOLATION interpolation;T_AVERAGING averaging;
        if(Z_min && Z_max) for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            TV X=iterator.Location()-dt*averaging.Face_To_Cell_Vector(grid,cell,face_velocities);
             Z(cell)=interpolation.Clamped_To_Array(grid,Z_ghost,X);
            VECTOR<T2,2> extrema=interpolation.Extrema_Clamped_To_Array_Cell(grid,*Z_min_ghost,*Z_max_ghost,X);
            (*Z_min)(cell)=extrema.x;(*Z_max)(cell)=extrema.y;}
        else for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            TV cell_velocity(averaging.Face_To_Cell_Vector(grid,cell,face_velocities));
            cell_velocity*=dt; // TODO: gcc 4.1.2 compiler bug workaround
            Z(cell)=interpolation.Clamped_To_Array(grid,Z_ghost,iterator.Location()-cell_velocity);}
    }
};

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION>
class ADVECTION_SEMI_LAGRANGIAN_TASK_FACE:public THREAD_QUEUE::TASK
{    
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;

public:

    const T_GRID& grid;
    T_FACE_ARRAYS_SCALAR& Z;
    const T_FACE_LOOKUP& Z_ghost;
    const T_FACE_LOOKUP& face_velocities;
    T_BOUNDARY& boundary;
    const T dt;
    const T time;
    const T_FACE_LOOKUP* Z_min_ghost;
    const T_FACE_LOOKUP* Z_max_ghost;
    T_FACE_ARRAYS_SCALAR* Z_min;
    T_FACE_ARRAYS_SCALAR* Z_max;
    RANGE<TV_INT> domain;
    int axis;

    ADVECTION_SEMI_LAGRANGIAN_TASK_FACE(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max,RANGE<TV_INT> domain,int axis)
        :grid(grid),Z(Z),Z_ghost(Z_ghost),face_velocities(face_velocities),boundary(boundary),dt(dt),time(time),Z_min_ghost(Z_min_ghost),Z_max_ghost(Z_max_ghost),Z_min(Z_min),Z_max(Z_max),domain(domain),axis(axis)
    {}
    
    void Run()
    {
        T_INTERPOLATION interpolation;T_AVERAGING averaging;
        if(Z_min && Z_max) for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
            TV_INT face=iterator.Face_Index();int axis=iterator.Axis();TV X=iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,axis,face,face_velocities);
            Z.Component(axis)(face)=interpolation.Clamped_To_Array_Face_Component(axis,grid,Z_ghost.Starting_Point_Face(axis,face),X);
            VECTOR<T,2> extrema=interpolation.Extrema_Clamped_To_Array_Face_Component(axis,grid,Z_min_ghost->Starting_Point_Face(axis,face),Z_max_ghost->Starting_Point_Face(axis,face),X);
            (*Z_min).Component(axis)(face)=extrema.x;(*Z_max).Component(axis)(face)=extrema.y;}
        else for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
            Z.Component(axis)(face)=interpolation.Clamped_To_Array_Face_Component(axis,grid,Z_ghost.Starting_Point_Face(axis,face),
                iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,axis,face,face_velocities));}
    }
};
}
#endif
#endif
