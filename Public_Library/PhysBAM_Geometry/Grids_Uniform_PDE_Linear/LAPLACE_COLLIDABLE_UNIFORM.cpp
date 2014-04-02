//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE_UNIFORM  
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_1D.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_2D.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LAPLACE_COLLIDABLE_UNIFORM<T_GRID>::
LAPLACE_COLLIDABLE_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input,THREAD_QUEUE* thread_queue)
    :BASE(grid_input,u_input,initialize_grid,enforce_compatibility_input,thread_queue),levelset_default(new T_LEVELSET(grid,phi_default))
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LAPLACE_COLLIDABLE_UNIFORM<T_GRID>::
LAPLACE_COLLIDABLE_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,T_LEVELSET& cell_centered_levelset,const bool initialize_grid,const bool multiphase_input,
    const bool enforce_compatibility_input,THREAD_QUEUE* thread_queue)
    :BASE(grid_input,u_input,initialize_grid,enforce_compatibility_input,thread_queue),levelset_default(new T_LEVELSET(grid,phi_default))
{
    levelset=&cell_centered_levelset;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LAPLACE_COLLIDABLE_UNIFORM<T_GRID>::
~LAPLACE_COLLIDABLE_UNIFORM()
{
    delete levelset_default;
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T_GRID> void LAPLACE_COLLIDABLE_UNIFORM<T_GRID>::
Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    BASE::Find_A_Part_Two(domain,A_array,b_array,cell_index_to_matrix_index);
    if(second_order_cut_cell_method) Apply_Second_Order_Cut_Cell_Method(domain,A_array,b_array,cell_index_to_matrix_index);
}
//#####################################################################
// Function Apply_Second_Order_Cut_Cell_Method
//#####################################################################
template<class T_GRID> void LAPLACE_COLLIDABLE_UNIFORM<T_GRID>::
Apply_Second_Order_Cut_Cell_Method(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    assert(levelset);
    TV minus_one_over_dx_squared=(T)-1*Inverse(grid.dX*grid.dX);
    for(int i=1;i<=TV::dimension;i++){
        RANGE<TV_INT> face_domain(domain);
        for(int axis=1;axis<=TV::dimension;axis++){
            if(face_domain.min_corner(axis)==grid.Domain_Indices(1).min_corner(axis)) face_domain.min_corner(axis)+=1;
            if(face_domain.max_corner(axis)==grid.Domain_Indices(1).max_corner(axis)) face_domain.max_corner(axis)-=1;}
        if(face_domain.min_corner(i)==grid.Domain_Indices().min_corner(i)) face_domain.min_corner(i)+=1;
        if(face_domain.max_corner(i)!=grid.Domain_Indices().max_corner(i)) face_domain.max_corner(i)+=1;
        for(FACE_ITERATOR iterator(grid,face_domain,i);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index(),cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
            if(!psi_N.Component(axis)(face) && LEVELSET_UTILITIES<T>::Interface(levelset->phi(cell1),levelset->phi(cell2))){
                T theta=LEVELSET_UTILITIES<T>::Theta(levelset->phi(cell1),levelset->phi(cell2));
                if(psi_D(cell1) && !psi_D(cell2) && domain.Lazy_Inside(cell2)){ // interface is to the negative side of second cell
                    int color=filled_region_colors(cell2);int matrix_index=cell_index_to_matrix_index(cell2);
                    T A_right_i=minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);b_array(color)(matrix_index)-=A_right_i*u(cell1);
                    A_right_i/=max((1-theta),second_order_cut_cell_threshold);
                    A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face);}  
                else if(psi_D(cell2) && !psi_D(cell1) && domain.Lazy_Inside(cell1)){ // interface is to the positive side of first cell
                    int color=filled_region_colors(cell1);int matrix_index=cell_index_to_matrix_index(cell1);
                    T A_right_i=minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);b_array(color)(matrix_index)-=A_right_i*u(cell2);
                    A_right_i/=max(theta,second_order_cut_cell_threshold);
                    A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face);}}}}
}
template<class T_GRID> void LAPLACE_COLLIDABLE_UNIFORM<T_GRID>::
Initialize_Grid(const T_GRID& mac_grid_input)
{
    assert(mac_grid_input.DX()==TV() || mac_grid_input.Is_MAC_Grid());
    BASE::Initialize_Grid(mac_grid_input);
    if(levelset == levelset_default) phi_default.Resize(grid.Domain_Indices(1));
    Set_Up_Second_Order_Cut_Cell_Method(second_order_cut_cell_method);
}
template<class T_GRID> void LAPLACE_COLLIDABLE_UNIFORM<T_GRID>::
Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input)
{
    second_order_cut_cell_method=use_second_order_cut_cell_method_input;
    if(second_order_cut_cell_method) u_interface.Resize(grid);
    else u_interface.Clean_Memory();
}
//#####################################################################
template class LAPLACE_COLLIDABLE_UNIFORM<GRID<VECTOR<float,1> > >;
template class LAPLACE_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> > >;
template class LAPLACE_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_COLLIDABLE_UNIFORM<GRID<VECTOR<double,1> > >;
template class LAPLACE_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> > >;
template class LAPLACE_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
}
