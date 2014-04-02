//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/POISSON_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
using namespace PhysBAM;
template<class T_GRID> POISSON_UNIFORM<T_GRID>::
POISSON_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input)
    :POISSON<T>(multiphase_input),LAPLACE_UNIFORM<T_GRID>(grid_input,u_input,false,enforce_compatibility_input)
{
    Initialize_Grid(grid_input);
}
template<class T_GRID> POISSON_UNIFORM<T_GRID>::
~POISSON_UNIFORM()
{
}
//#####################################################################
// Function Find_Variable_beta
//#####################################################################
// only set up for Dirichlet boundary conditions - doesn't work for jump conditons yet 
template<class T_GRID> void POISSON_UNIFORM<T_GRID>::
Find_Variable_beta()
{
    if(beta_given_on_faces) return; // beta_right, beta_top, beta_back already set
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) 
        beta_face.Component(iterator.Axis())(iterator.Face_Index())=(T).5*(variable_beta(iterator.First_Cell_Index())+variable_beta(iterator.Second_Cell_Index())); 
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T_GRID> void POISSON_UNIFORM<T_GRID>::
Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    TV one_over_dx2=Inverse(grid.dX*grid.dX);
    TV_INT grid_counts=grid.counts;
    if(use_weighted_divergence)
        for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
            int color=filled_region_colors(iterator.Cell_Index());
            if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){const TV_INT& cell_index=iterator.Cell_Index();
                int matrix_index=cell_index_to_matrix_index(cell_index);
                SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(filled_region_colors(cell_index));VECTOR_ND<T>& b=b_array(filled_region_colors(cell_index));b(matrix_index)=f(cell_index);
                T diagonal=0;
                for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
                    if(filled_region_colors.Valid_Index(cell_index-offset)){
                        if(!psi_N.Component(axis)(cell_index)){
                            T element=divergence_face_weights.Component(axis)(cell_index)*beta_face.Component(axis)(cell_index)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index-offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index-offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index-offset)) b(matrix_index)-=element*u(cell_index-offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),element);}}
                    if(filled_region_colors.Valid_Index(cell_index+offset)){
                        if(!psi_N.Component(axis)(cell_index+offset)){
                            T element=divergence_face_weights.Component(axis)(cell_index+offset)*beta_face.Component(axis)(cell_index+offset)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index+offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index+offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index+offset)) b(matrix_index)-=element*u(cell_index+offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),element);}}}
                A.Set_Element(matrix_index,matrix_index,diagonal);}} // set diagonal and right hand side
    else
        for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
            int color=filled_region_colors(iterator.Cell_Index());
            if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){const TV_INT& cell_index=iterator.Cell_Index();
                int matrix_index=cell_index_to_matrix_index(cell_index);
                SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(filled_region_colors(cell_index));VECTOR_ND<T>& b=b_array(filled_region_colors(cell_index));b(matrix_index)=f(cell_index);
                T diagonal=0;
                for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
                    if(filled_region_colors.Valid_Index(cell_index-offset)){
                        if(!psi_N.Component(axis)(cell_index)){T element=beta_face.Component(axis)(cell_index)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index-offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index-offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index-offset)) b(matrix_index)-=element*u(cell_index-offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),element);}}
                    if(filled_region_colors.Valid_Index(cell_index+offset)){
                        if(!psi_N.Component(axis)(cell_index+offset)){T element=beta_face.Component(axis)(cell_index+offset)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index+offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index+offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index+offset)) b(matrix_index)-=element*u(cell_index+offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),element);}}}
                A.Set_Element(matrix_index,matrix_index,diagonal);}} // set diagonal and right hand side
}
template<class T_GRID> void POISSON_UNIFORM<T_GRID>::
Initialize_Grid(const T_GRID& grid_input)
{
    LAPLACE_UNIFORM<T_GRID>::Initialize_Grid(grid_input); // laplace will call Set_Up_Second_Order_Cut_Cell_Method to resize those arrays
    beta_face.Resize(grid,1); // TODO: check if can get rid of extra ghost layer
    if(use_variable_beta)variable_beta.Resize(grid.Domain_Indices(1));
}
//#####################################################################
template class POISSON_UNIFORM<GRID<VECTOR<float,1> > >;
template class POISSON_UNIFORM<GRID<VECTOR<float,2> > >;
template class POISSON_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POISSON_UNIFORM<GRID<VECTOR<double,1> > >;
template class POISSON_UNIFORM<GRID<VECTOR<double,2> > >;
template class POISSON_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
