//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
template<class T,class T_GRID> typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_SCALAR LEVELSET<T,T_GRID>::interpolation_default;
template<class T,class T_GRID> typename REBIND<typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_SCALAR,typename T_GRID::VECTOR_T>::TYPE LEVELSET<T,T_GRID>::normal_interpolation_default;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LEVELSET_UNIFORM<T_GRID>::
LEVELSET_UNIFORM(T_GRID& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells_input)
    :grid(grid_input),phi(phi_input),normals(0),curvature(0),cell_range(0),thread_queue(0),number_of_ghost_cells(number_of_ghost_cells_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LEVELSET_UNIFORM<T_GRID>::
~LEVELSET_UNIFORM()
{
    delete normals;delete curvature;delete cell_range;
}
//#####################################################################
// Function Collision_Aware_Phi
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_UNIFORM<T_GRID>::
Collision_Aware_Phi(const TV& location) const
{
    assert(collision_body_list);
    return LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value).Clamped_To_Array(grid,phi,location);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_UNIFORM<T_GRID>::
CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const
{
    T dt_convection=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        T local_V_norm=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++)
            local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}
    T dt_curvature=(curvature_motion && T_GRID::dimension>1)?sigma*2*grid.one_over_dX.Magnitude_Squared():0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_UNIFORM<T_GRID>::
CFL(const T_ARRAYS_VECTOR& velocity) const
{
    T dt_convection=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        dt_convection=max(dt_convection,TV::Dot_Product(velocity(cell),grid.one_over_dX));}
    T dt_curvature=(curvature_motion && T_GRID::dimension>1)?sigma*2*grid.one_over_dX.Magnitude_Squared():0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Iterative_Find_Interface
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_UNIFORM<T_GRID>::
Iterative_Find_Interface(TV left,TV right,const int iterations) const
{
    T phi_left=Phi(left),phi_right=Phi(right);
    TV interface=LINEAR_INTERPOLATION<T,TV>::Linear(left,right,LEVELSET_UTILITIES<T>::Theta(phi_left,phi_right));
    int phi_left_sign=(phi_left<=0?-1:1),phi_right_sign=(phi_right<=0?-1:1);
    for(int i=1;i<=iterations;i++){
        T phi=Phi(interface);
        int phi_sign=(phi<=0?-1:1);
        if(phi_left_sign*phi_sign<0){
            right=interface;phi_right=phi;phi_right_sign=phi_sign;
            interface=LINEAR_INTERPOLATION<T,TV>::Linear(left,interface,LEVELSET_UTILITIES<T>::Theta(phi_left,phi));}
        else if(phi_right_sign*phi_sign<0){
            left=interface;phi_left=phi;phi_left_sign=phi_sign;
            interface=LINEAR_INTERPOLATION<T,TV>::Linear(interface,right,LEVELSET_UTILITIES<T>::Theta(phi,phi_right));}
        else break;}
    return interface;
}
//#####################################################################
// Function Compute_Gradient
//#####################################################################
template<class T_GRID> void LEVELSET_UNIFORM<T_GRID>::
Compute_Gradient(T_ARRAYS_VECTOR& gradient,const T time) const
{
    TV one_over_two_dx=(T).5*grid.one_over_dX;
    int ghost_cells=3;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    gradient.Resize(grid.Domain_Indices(ghost_cells-1));
    for(CELL_ITERATOR iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        for(int axis=1;axis<=TV::dimension;axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            gradient(cell_index)(axis)=(phi_ghost(cell_index+axis_vector)-phi_ghost(cell_index-axis_vector))*one_over_two_dx(axis);}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class LEVELSET_UNIFORM<GRID<VECTOR<T,d> > >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
