//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_PDE_Linear/POISSON_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_NXN.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_beta_And_Add_Jumps_To_b
//#####################################################################
template<class T_GRID> void POISSON_DYADIC<T_GRID>::
Compute_beta_And_Add_Jumps_To_b(const T dt,const T time)
{
}
//#####################################################################
// Function Find_Constant_beta
//#####################################################################
// only set up for jump conditons - doesn't work for Dirichlet boundary conditions yet
template<class T_GRID> void POISSON_DYADIC<T_GRID>::
Find_Constant_beta()
{      
    PHYSBAM_FATAL_ERROR("Not Implemented");
}
//#####################################################################
// Function Find_Variable_beta
//#####################################################################
// only set up for Dirichlet boundary conditions - doesn't work for jump conditons yet
template<class T_GRID> void POISSON_DYADIC<T_GRID>::
Find_Variable_beta()
{
    if(beta_given_on_faces) return; // beta_face already set
    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
        beta_face(iterator.Face_Index())=(T).5*(variable_beta(iterator.First_Cell_Index())+variable_beta(iterator.Second_Cell_Index()));
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T_GRID> void POISSON_DYADIC<T_GRID>::
Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index)
{
    ARRAY<ARRAY<int> > row_counts(A_array.m,false);
    ARRAY<ARRAY<T> > row_sum(A_array.m,false);
    for(int i=1;i<=A_array.m;i++){
        row_counts(i).Resize(filled_region_cell_count(i));
        row_sum(i).Resize(filled_region_cell_count(i));
        b_array(i).Resize(filled_region_cell_count(i));}
    Find_A_Part_One(filled_region_cell_count,cell_index_to_matrix_index,row_counts,row_sum);
    for(int i=1;i<=A_array.m;i++) A_array(i).Set_Row_Lengths(row_counts(i));
    Find_A_Part_Two(A_array,b_array,cell_index_to_matrix_index,row_sum);
}
//#####################################################################
// Function Find_A_Part_One
//#####################################################################
template<class T_GRID> void POISSON_DYADIC<T_GRID>::
Find_A_Part_One(const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& row_counts_array,ARRAY<ARRAY<T> >& row_sum_array)
{
    // TODO(jontg): periodic boundary conditions
    // TODO(jontg): touches dirichlet region / solve_neumann_regions
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){INDEX cell_index=iterator.Cell_Index();
        int color=filled_region_colors(cell_index);
        if(!psi_D(cell_index) && color!=-1) row_counts_array(color)(cell_index_to_matrix_index(cell_index))=1;}

    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        if(psi_N(iterator.Face_Index())) continue;
        INDEX left_cell_index=iterator.First_Cell_Index(),right_cell_index=iterator.Second_Cell_Index(); int axis=iterator.Axis();
        int left_color=filled_region_colors(left_cell_index),right_color=filled_region_colors(right_cell_index),active_color=left_color;
        if(active_color==-1 && right_color==-1) continue; // both uncolorable
        else if(active_color != -1 && right_color != -1 && active_color != right_color) continue; // different colors
        else if(active_color==-1) active_color=right_color;

        ARRAY<T>& row_sum=row_sum_array(active_color);
        ARRAY<int>& row_counts=row_counts_array(active_color);

        const CELL* smaller_cell=iterator.Deepest_Cell();const CELL* larger_cell=iterator.Other_Cell();
        T smaller_dx=smaller_cell->DX()[axis+1],larger_dx=larger_cell->DX()[axis+1];
        int smaller_index=cell_index_to_matrix_index(smaller_cell->Cell()),larger_index=cell_index_to_matrix_index(larger_cell->Cell());

        T distance=(T).5*(smaller_dx+larger_dx),entry=beta_face(iterator.Face_Index())*smaller_cell->Face_Size()/distance;
        if(!psi_D(smaller_cell->Cell()) && !psi_D(larger_cell->Cell())){ // if neither cell is dirichlet
            row_sum(larger_index)+=entry;++row_counts(larger_index);
            row_sum(smaller_index)+=entry;++row_counts(smaller_index);}
        else if(!psi_D(larger_cell->Cell())) row_sum(larger_index)+=entry;
        else if(!psi_D(smaller_cell->Cell())) row_sum(smaller_index)+=entry;}
}
//#####################################################################
// Function Find_A_Part_Two
//#####################################################################
template<class T_GRID> void POISSON_DYADIC<T_GRID>::
Find_A_Part_Two(ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<T> >& row_sum_array)
{
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){int cell_index=iterator.Cell_Index();int matrix_index=cell_index_to_matrix_index(cell_index);
        int color=filled_region_colors(cell_index);
        if(psi_D(cell_index) || color==-1) continue;
        A_array(color).Set_Element(matrix_index,matrix_index,-row_sum_array(color)(matrix_index));
        b_array(color)(matrix_index)+=iterator.Cell_Pointer()->Cell_Size()*f(cell_index);}

    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        if(psi_N(iterator.Face_Index())) continue;
        INDEX left_cell_index=iterator.First_Cell_Index(),right_cell_index=iterator.Second_Cell_Index(); int axis=iterator.Axis();
        int left_color=filled_region_colors(left_cell_index),right_color=filled_region_colors(right_cell_index),active_color=left_color;
        if(active_color==-1 && right_color==-1) continue; // both uncolorable
        else if(active_color != -1 && right_color != -1 && active_color != right_color) continue; // different colors
        else if(active_color==-1) active_color=right_color;

        const CELL* smaller_cell=iterator.Deepest_Cell();
        const CELL* larger_cell=iterator.Other_Cell();
        T smaller_dx=smaller_cell->DX()[axis+1],larger_dx=larger_cell->DX()[axis+1];
        int smaller_index=cell_index_to_matrix_index(smaller_cell->Cell()),larger_index=cell_index_to_matrix_index(larger_cell->Cell());

        T distance=(T).5*(smaller_dx+larger_dx),entry=beta_face(iterator.Face_Index())*smaller_cell->Face_Size()/distance;
        if(!psi_D(smaller_cell->Cell()) && !psi_D(larger_cell->Cell())) // if neither cell is dirichlet
            A_array(active_color).Add_Symmetric_Elements(smaller_index,larger_index,entry);
        else if(!psi_D(larger_cell->Cell())) b_array(active_color)(larger_index)-=entry*u(smaller_cell->Cell());
        else if(!psi_D(smaller_cell->Cell())) b_array(active_color)(smaller_index)-=entry*u(larger_cell->Cell());}
}
//#####################################################################
#if COMPILE_WITH_BINTREE_SUPPORT
template class POISSON_DYADIC<BINTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POISSON_DYADIC<BINTREE_GRID<double> >;
#endif
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class POISSON_DYADIC<OCTREE_GRID<float> >;
template class POISSON_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POISSON_DYADIC<OCTREE_GRID<double> >;
template class POISSON_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif

#endif
