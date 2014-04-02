//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_MACCORMACK_UNIFORM
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2,class T_NESTED_ADVECTION> ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
ADVECTION_MACCORMACK_UNIFORM(T_NESTED_ADVECTION& nested_advection_input,const T_ARRAYS_BOOL* node_mask_input,const T_ARRAYS_BOOL* cell_mask_input,const T_FACE_ARRAYS_BOOL* face_mask_input,THREAD_QUEUE* thread_queue_input)
    :node_mask(node_mask_input),cell_mask(cell_mask_input),face_mask(face_mask_input),nested_advection(nested_advection_input),clamp_extrema(true),ensure_second_order(false),thread_queue(thread_queue_input)
{
}
//#####################################################################
// Function Update_Advection_Equation_Node
//#####################################################################
template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost_input,const T_ARRAYS_T2* Z_max_ghost_input,T_ARRAYS_T2* Z_min_input,T_ARRAYS_T2* Z_max_input)
{
    assert(node_mask && !Z_min_input && !Z_max_input);
    T_ARRAYS_VECTOR negative_V(grid.Domain_Indices(0),false);T_ARRAYS_VECTOR::Copy((T)-1,V,negative_V);
    T_ARRAYS_T2 Z_temp=Z,Z_min_ghost=Z_ghost,Z_max_ghost=Z_ghost,Z_min(grid.Domain_Indices(),false),Z_max(grid.Domain_Indices(),false);
    int number_of_ghost_cells=Z.Domain_Indices().Minimum_Corner()(1)-Z_ghost.Domain_Indices().Minimum_Corner()(1);
    nested_advection.Update_Advection_Equation_Node(grid,Z_temp,Z_ghost,V,boundary,dt,time,&Z_min_ghost,&Z_max_ghost,&Z_min,&Z_max); // Z_temp now has time n+1 data
    T_ARRAYS_T2 Z_forward_ghost(grid.Domain_Indices(number_of_ghost_cells),false);boundary.Fill_Ghost_Cells(grid,Z_temp,Z_forward_ghost,dt,time+dt,number_of_ghost_cells);
    boundary.Fill_Ghost_Cells(grid,Z_min,Z_min_ghost,dt,time+dt,number_of_ghost_cells);boundary.Fill_Ghost_Cells(grid,Z_max,Z_max_ghost,dt,time+dt,number_of_ghost_cells);
    nested_advection.Update_Advection_Equation_Node(grid,Z_temp,Z_forward_ghost,negative_V,boundary,dt,time+dt,&Z_min_ghost,&Z_max_ghost,&Z_min,&Z_max); // Z_temp has time n data
    RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::All_Ones_Vector();
    DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>,typename T_GRID::VECTOR_T> thread_iterator(domain,thread_queue);
    if(clamp_extrema)
        thread_iterator.template Run<const T_GRID&,T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&>(*this,&ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::Apply_Clamped_Extrema_Limiter_Node_Threaded,grid,Z,Z_forward_ghost,Z_temp,Z_min,Z_max);
    else
        thread_iterator.template Run<const T_GRID&,T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&>(*this,&ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::Apply_Reversion_Limiter_Node_Threaded,grid,Z,Z_forward_ghost,Z_temp,Z_min,Z_max);
}
//#####################################################################
// Function Update_Advection_Equation_Cell_Lookup
//#####################################################################
template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const FACE_LOOKUP_UNIFORM<T_GRID>& face_velocities,T_BOUNDARY_T2& boundary,
    const T dt,const T time,const T_ARRAYS_T2* Z_min_ghost_input,const T_ARRAYS_T2* Z_max_ghost_input,T_ARRAYS_T2* Z_min_input,T_ARRAYS_T2* Z_max_input)
{
    assert(cell_mask && !Z_min_input && !Z_max_input);
    T_FACE_ARRAYS_SCALAR negative_V(face_velocities.V_face.Domain_Indices());T_FACE_ARRAYS_SCALAR::Copy(-1,face_velocities.V_face,negative_V);
    T_ARRAYS_T2 Z_temp=Z,Z_min_ghost=Z_ghost,Z_max_ghost=Z_ghost,Z_min(grid.Domain_Indices(),false),Z_max(grid.Domain_Indices(),false);
    int number_of_ghost_cells=Z.Domain_Indices().Minimum_Corner()(1)-Z_ghost.Domain_Indices().Minimum_Corner()(1);
    nested_advection.Update_Advection_Equation_Cell(grid,Z_temp,Z_ghost,face_velocities.V_face,boundary,dt,time,&Z_min_ghost,&Z_max_ghost,&Z_min,&Z_max); // Z_temp now has time n+1 data
    T_ARRAYS_T2 Z_forward_ghost(grid.Domain_Indices(number_of_ghost_cells),false);boundary.Fill_Ghost_Cells(grid,Z_temp,Z_forward_ghost,dt,time+dt,number_of_ghost_cells);
    boundary.Fill_Ghost_Cells(grid,Z_min,Z_min_ghost,dt,time+dt,number_of_ghost_cells);boundary.Fill_Ghost_Cells(grid,Z_max,Z_max_ghost,dt,time+dt,number_of_ghost_cells);
    nested_advection.Update_Advection_Equation_Cell(grid,Z_temp,Z_forward_ghost,negative_V,boundary,dt,time+dt,&Z_min_ghost,&Z_max_ghost,&Z_min,&Z_max); // Z_temp has time n data
    DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>,typename T_GRID::VECTOR_T> thread_iterator(grid.Domain_Indices(),thread_queue);
    if(clamp_extrema)
        thread_iterator.template Run<const T_GRID&,T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&>(*this,&ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::Apply_Clamped_Extrema_Limiter_Cell_Threaded,grid,Z,Z_forward_ghost,Z_temp,Z_min,Z_max);
    else
        thread_iterator.template Run<const T_GRID&,T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&,const T_ARRAYS_T2&>(*this,&ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::Apply_Reversion_Limiter_Cell_Threaded,grid,Z,Z_forward_ghost,Z_temp,Z_min,Z_max);
}
//#####################################################################
// Function Update_Advection_Equation_Face_Lookup
//#####################################################################
template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const FACE_LOOKUP_UNIFORM<T_GRID>& Z_ghost,const FACE_LOOKUP_UNIFORM<T_GRID>& face_velocities,
    T_BOUNDARY& boundary,const T dt,const T time,const FACE_LOOKUP_UNIFORM<T_GRID>* Z_min_ghost_input,const FACE_LOOKUP_UNIFORM<T_GRID>* Z_max_ghost_input,
    T_FACE_ARRAYS_SCALAR* Z_min_input,T_FACE_ARRAYS_SCALAR* Z_max_input)
{
    assert(!ensure_second_order || !clamp_extrema);
    assert(face_mask && !Z_min_input && !Z_max_input);
    T_FACE_ARRAYS_SCALAR Z_temp=Z;
    int number_of_ghost_cells=Z_ghost.Number_Of_Ghost_Cells();
    if(!ensure_second_order){
        T_FACE_ARRAYS_SCALAR Z_min_ghost=Z_ghost.V_face,Z_max_ghost=Z_ghost.V_face,Z_min(grid,0,false),Z_max(grid,0,false);
        nested_advection.Update_Advection_Equation_Face(grid,Z_temp,Z_ghost.V_face,face_velocities.V_face,boundary,dt,time,&Z_min_ghost,&Z_max_ghost,&Z_min,&Z_max); // Z_temp now has time n+1 data
        T_FACE_ARRAYS_SCALAR Z_forward_ghost(grid,number_of_ghost_cells,false);boundary.Fill_Ghost_Cells_Face(grid,Z_temp,Z_forward_ghost,time+dt,number_of_ghost_cells);
        boundary.Fill_Ghost_Cells_Face(grid,Z_min,Z_min_ghost,time+dt,number_of_ghost_cells);boundary.Fill_Ghost_Cells_Face(grid,Z_max,Z_max_ghost,time+dt,number_of_ghost_cells);
        nested_advection.Update_Advection_Equation_Face(grid,Z_temp,Z_forward_ghost,face_velocities.V_face,boundary,-dt,time+dt,&Z_min_ghost,&Z_max_ghost,&Z_min,&Z_max); // Z_temp has time n data
        if(clamp_extrema)
            for(int axis=1;axis<=T_GRID::dimension;++axis){
                RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)++;
                DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>,typename T_GRID::VECTOR_T> thread_iterator(domain,thread_queue);
                thread_iterator.template Run<const T_GRID&,int,T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&>(*this,&ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::Apply_Clamped_Extrema_Limiter_Face_Threaded,grid,axis,Z,Z_forward_ghost,Z_temp,Z_min,Z_max);}
        else
            for(int axis=1;axis<=T_GRID::dimension;++axis){
                RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)++;
                DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>,typename T_GRID::VECTOR_T> thread_iterator(domain,thread_queue);
                thread_iterator.template Run<const T_GRID&,int,T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&>(*this,&ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::Apply_Reversion_Limiter_Face_Threaded,grid,axis,Z,Z_forward_ghost,Z_temp,Z_min,Z_max);}}
    else{
        nested_advection.Update_Advection_Equation_Face(grid,Z_temp,Z_ghost.V_face,face_velocities.V_face,boundary,dt,time,0,0,0,0); // Z_temp now has time n+1 data
        T_FACE_ARRAYS_SCALAR Z_forward_ghost(grid,number_of_ghost_cells,false);boundary.Fill_Ghost_Cells_Face(grid,Z_temp,Z_forward_ghost,time+dt,number_of_ghost_cells);
        nested_advection.Update_Advection_Equation_Face(grid,Z_temp,Z_forward_ghost,face_velocities.V_face,boundary,-dt,time+dt,0,0,0,0); // Z_temp has time n data
        for(int axis=1;axis<=T_GRID::dimension;++axis){
            RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)++;
            DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>,typename T_GRID::VECTOR_T> thread_iterator(domain,thread_queue);
            thread_iterator.template Run<const T_GRID&,int,T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&,const T_FACE_ARRAYS_SCALAR&>(*this,&ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::Apply_Second_Order_Update_Face_Threaded,grid,axis,Z,Z_forward_ghost,Z_temp);}}
}
template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Apply_Clamped_Extrema_Limiter_Node_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_forward_ghost,const T_ARRAYS_T2& Z_backward_ghost,const T_ARRAYS_T2& Z_min,const T_ARRAYS_T2& Z_max)
{
    for(NODE_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Node_Index();
        if((*node_mask)(index)) Z(index)=clamp(Z_forward_ghost(index)+(T).5*(Z(index)-Z_backward_ghost(index)),Z_min(index),Z_max(index));
        else Z(index)=Z_forward_ghost(index);}
}

template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Apply_Clamped_Extrema_Limiter_Cell_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z, const T_ARRAYS_T2& Z_forward_ghost, const T_ARRAYS_T2& Z_backward_ghost, const T_ARRAYS_T2& Z_min, const T_ARRAYS_T2& Z_max)
{
    for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Cell_Index();
        if((*cell_mask)(index)) Z(index)=clamp(Z_forward_ghost(index)+(T).5*(Z(index)-Z_backward_ghost(index)),Z_min(index),Z_max(index));
        else Z(index)=Z_forward_ghost(index);}
}

template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Apply_Clamped_Extrema_Limiter_Face_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost,const T_FACE_ARRAYS_SCALAR& Z_min,const T_FACE_ARRAYS_SCALAR& Z_max)
{
    for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Face_Index();int axis=iterator.Axis();
        if((*face_mask)(axis,index)) Z(axis,index)=clamp(Z_forward_ghost(axis,index)+(T).5*(Z(axis,index)-Z_backward_ghost(axis,index)),Z_min(axis,index),Z_max(axis,index));
        else Z(axis,index)=Z_forward_ghost(axis,index);}
}

template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Apply_Reversion_Limiter_Node_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_forward_ghost,const T_ARRAYS_T2& Z_backward_ghost,const T_ARRAYS_T2& Z_min,const T_ARRAYS_T2& Z_max)
{
    for(NODE_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Node_Index();
        if((*node_mask)(index)){
            Z(index)=Z_forward_ghost(index)+(T).5*(Z(index)-Z_backward_ghost(index));
            if(!in_bounds(Z(index),Z_min(index),Z_max(index))) Z(index)=Z_forward_ghost(index);}
        else Z(index)=Z_forward_ghost(index);}
}

template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Apply_Reversion_Limiter_Cell_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_forward_ghost,const T_ARRAYS_T2& Z_backward_ghost,const T_ARRAYS_T2& Z_min,const T_ARRAYS_T2& Z_max)
{
    for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Cell_Index();
        if((*cell_mask)(index)){
            Z(index)=Z_forward_ghost(index)+(T).5*(Z(index)-Z_backward_ghost(index));
            if(!in_bounds(Z(index),Z_min(index),Z_max(index))) Z(index)=Z_forward_ghost(index);}
        else Z(index)=Z_forward_ghost(index);}
}

template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Apply_Reversion_Limiter_Face_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost,const T_FACE_ARRAYS_SCALAR& Z_min,const T_FACE_ARRAYS_SCALAR& Z_max)
{
    for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Face_Index();int axis=iterator.Axis();
        if((*face_mask)(axis,index)){
            Z(axis,index)=Z_forward_ghost(axis,index)+(T).5*(Z(axis,index)-Z_backward_ghost(axis,index));
            if(!in_bounds(Z(axis,index),Z_min(axis,index),Z_max(axis,index))) Z(axis,index)=Z_forward_ghost(axis,index);}
        else Z(axis,index)=Z_forward_ghost(axis,index);}
}

template<class T_GRID,class T2,class T_NESTED_ADVECTION> void ADVECTION_MACCORMACK_UNIFORM<T_GRID,T2,T_NESTED_ADVECTION>::
Apply_Second_Order_Update_Face_Threaded(RANGE<TV_INT>& domain,const T_GRID& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost)
{
    for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> index=iterator.Full_Index();
        if((*face_mask)(index)) Z(index)=Z_forward_ghost(index)+(T).5*(Z(index)-Z_backward_ghost(index));
        else Z(index)=Z_forward_ghost(index);}
}

//#####################################################################
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,1> >,float,ADVECTION<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,1> >,float,ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,1> >,float,ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,2> >,float,ADVECTION<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,2> >,float,ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,2> >,float,ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,3> >,float,ADVECTION<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,3> >,float,ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<VECTOR<float,3> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<float,3> >,float,ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<VECTOR<float,3> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,1> >,double,ADVECTION<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,1> >,double,ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,1> >,double,ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,2> >,double,ADVECTION<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,2> >,double,ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,2> >,double,ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,3> >,double,ADVECTION<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,3> >,double,ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<VECTOR<double,3> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > > >;
template class ADVECTION_MACCORMACK_UNIFORM<GRID<VECTOR<double,3> >,double,ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<VECTOR<double,3> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > > >;
#endif
