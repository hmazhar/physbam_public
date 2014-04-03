//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
//#####################################################################
// Function Compute_Level_Set
//#####################################################################
template<class T> void LEVELSET_MAKER_UNIFORM_2D<T>::
Compute_Level_Set(SEGMENTED_CURVE_2D<T>& curve,GRID<TV>& grid,int ghost_cells,ARRAY<T,TV_INT>& phi)
{
    phi.Fill(FLT_MAX);
    T dx=grid.dX.Max();

    ARRAY<bool,TV_INT> done(grid.Domain_Indices(ghost_cells+1));
    for(int i=1;i<=curve.mesh.elements.m;i++){
        SEGMENT_2D<T> segment(curve.particles.X(curve.mesh.elements(i).x),curve.particles.X(curve.mesh.elements(i).y));
        RANGE<TV_INT> box(grid.Cell(segment.x1,3));
        box.Enlarge_To_Include_Point(grid.Cell(segment.x2,3));
        box=box.Intersect(box,grid.Domain_Indices(ghost_cells-1));
        if(box.Empty()) continue;
        for(UNIFORM_ARRAY_ITERATOR<TV::m> it(box.Thickened(1));it.Valid();it.Next()){
            TV X=grid.X(it.index);
            T dist=segment.Distance_From_Point_To_Segment(X);
            if(dist<abs(phi(it.index))+dx*1e-4 && dist<dx){
                bool new_sign=TV::Dot_Product(X-segment.x1,segment.Normal())<0;
                if(abs(dist-abs(phi(it.index)))<dx*1e-4 && new_sign != (phi(it.index)<0))
                    new_sign=curve.Inside(grid.X(it.index));
                if(abs(dist)<abs(phi(it.index))) phi(it.index)=dist;
                phi(it.index)=abs(phi(it.index))*(new_sign?-1:1);
                done(it.index)=true;}}}

    ARRAY<TV_INT> todo,next_todo;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        if(phi(it.index)!=FLT_MAX)
            todo.Append(it.index);
    for(int layer=1;todo.m;layer++){
        while(todo.m){
            TV_INT index=todo.Pop();
            T next=sign(phi(index))*layer*dx;
            Compute_Level_Set_Helper(index+TV_INT(1,0),next,next_todo,phi);
            Compute_Level_Set_Helper(index-TV_INT(1,0),next,next_todo,phi);
            Compute_Level_Set_Helper(index+TV_INT(0,1),next,next_todo,phi);
            Compute_Level_Set_Helper(index-TV_INT(0,1),next,next_todo,phi);}
        todo.Exchange(next_todo);}
/*
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        VECTOR<T,3> color=phi(it.index)>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0);
        if(done(it.index)) color/=(T)3;
        Add_Debug_Particle(grid.X(it.index),color);}*/
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Compute Level Set",0,1);

    LEVELSET_2D<GRID<TV> > levelset(grid,phi);
    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(levelset,ghost_cells);
    fmm.Fast_Marching_Method(phi,done);
}
//#####################################################################
// Function Compute_Level_Set_Helper
//#####################################################################
template<class T> void LEVELSET_MAKER_UNIFORM_2D<T>::
Compute_Level_Set_Helper(const TV_INT& index,T next,ARRAY<TV_INT>& next_todo,ARRAY<T,TV_INT>& phi)
{
    if(!phi.Valid_Index(index)) return;
    T& p=phi(index);
    if(p!=FLT_MAX) return;
    p=next;
    next_todo.Append(index);
}
template class LEVELSET_MAKER_UNIFORM_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_MAKER_UNIFORM_2D<double>;
#endif
