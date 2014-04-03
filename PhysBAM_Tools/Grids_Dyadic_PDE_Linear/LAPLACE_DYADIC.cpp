//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/FLOOD_FILL_BINTREE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/FLOOD_FILL_OCTREE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/FLOOD_FILL_QUADTREE.h>
#include <PhysBAM_Tools/Grids_Dyadic_PDE_Linear/LAPLACE_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_NXN.h>
using namespace PhysBAM;
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
Solve(const T time,const bool solution_regions_already_computed)
{
    if(!solution_regions_already_computed) Find_Solution_Regions();

    // build all the matrices
    ARRAY<ARRAY<int> > matrix_index_to_cell_index_array(number_of_regions);ARRAY<int> cell_index_to_matrix_index(grid.number_of_cells);
    ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-1,number_of_regions);
    // for(int i=1;i<=grid.number_of_cells;i++) 
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,1);iterator.Valid();iterator.Next()) filled_region_cell_count(filled_region_colors(iterator.Cell_Index()))++;
    ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array(number_of_regions);ARRAY<VECTOR_ND<T> > b_array(number_of_regions);
    for(int color=1;color<=number_of_regions;color++) if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));}
    filled_region_cell_count.Fill(0);
    // for(int i=1;i<=grid.number_of_cells;i++){
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,1);iterator.Valid();iterator.Next()){
        int i=iterator.Cell_Index();
        int color=filled_region_colors(i);
        if(color>0&&(filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            cell_index_to_matrix_index(i)=++filled_region_cell_count(color);
            matrix_index_to_cell_index_array(color)(filled_region_cell_count(color))=i;}
    }
    RANGE<TV_INT> dummy_domain(TV_INT(-1),TV_INT(1));
    Find_A(dummy_domain,A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);

    for(int color=1;color<=number_of_regions;color++) if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        if(enforce_compatibility) pcg.Enforce_Compatibility(!filled_region_touches_dirichlet(color));
        Solve_Subregion(matrix_index_to_cell_index_array(color),A_array(color),b_array(color));}
    if(!solve_neumann_regions) //for(int i=1;i<=grid.number_of_cells;i++){
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,0);iterator.Valid();iterator.Next()){int i=iterator.Cell_Index();
        int filled_region_color=filled_region_colors(i);if(filled_region_color>0 && !filled_region_touches_dirichlet(filled_region_color)) u(i)=0;}
}
//#####################################################################
// Function Solve_Subregion
//#####################################################################
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
Solve_Subregion(ARRAY<int>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b)
{
    int number_of_unknowns=matrix_index_to_cell_index.m;
    A.Negate();b*=(T)-1;
    VECTOR_ND<T> x(number_of_unknowns),q,s,r,k,z;
    for(int i=1;i<=number_of_unknowns;i++) x(i)=u(matrix_index_to_cell_index(i));
    Find_Tolerance(b); // needs to happen after b is completely set up
    if(pcg.show_results) LOG::cout << "solving " << number_of_unknowns << " cells to tolerance " << tolerance << std::endl;
    pcg.Solve(A,x,b,q,s,r,k,z,tolerance);
    for(int i=1;i<=number_of_unknowns;i++)u(matrix_index_to_cell_index(i))=x(i);
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index)
{
    ARRAY<ARRAY<int> > row_counts(A_array.m,false);
    ARRAY<ARRAY<T> > row_sum(A_array.m,false);
    for(int i=1;i<=A_array.m;i++){
        row_counts(i).Resize(filled_region_cell_count(i),false,false);
        row_sum(i).Resize(filled_region_cell_count(i),false,false);
        b_array(i).Resize(filled_region_cell_count(i));}
    Find_A_Part_One(filled_region_cell_count,cell_index_to_matrix_index,row_counts,row_sum);
    for(int i=1;i<=A_array.m;i++) A_array(i).Set_Row_Lengths(row_counts(i));
    Find_A_Part_Two(A_array,b_array,cell_index_to_matrix_index,row_sum);
}
//#####################################################################
// Function Find_A_Part_One
//#####################################################################
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
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

        T distance=(T).5*(smaller_dx+larger_dx),entry=smaller_cell->Face_Size()/distance;
        if(!psi_D(smaller_cell->Cell()) && !psi_D(larger_cell->Cell())){ // if neither cell is dirichlet
            row_sum(larger_index)+=entry;++row_counts(larger_index);
            row_sum(smaller_index)+=entry;++row_counts(smaller_index);}
        else if(!psi_D(larger_cell->Cell())) row_sum(larger_index)+=entry;
        else if(!psi_D(smaller_cell->Cell())) row_sum(smaller_index)+=entry;}
}
//#####################################################################
// Function Find_A_Part_Two
//#####################################################################
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
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

        T distance=(T).5*(smaller_dx+larger_dx),entry=smaller_cell->Face_Size()/distance;
        if(!psi_D(smaller_cell->Cell()) && !psi_D(larger_cell->Cell())) // if neither cell is dirichlet
            A_array(active_color).Add_Symmetric_Elements(smaller_index,larger_index,entry);
        else if(!psi_D(larger_cell->Cell())) b_array(active_color)(larger_index)-=entry*u(smaller_cell->Cell());
        else if(!psi_D(smaller_cell->Cell())) b_array(active_color)(smaller_index)-=entry*u(larger_cell->Cell());}
}
//#####################################################################
// Function Find_Solution_Regions
//#####################################################################
template<class T_GRID> static void Mark_Interior_Cells(void* data,const typename T_GRID::CELL* cell)
{
    LAPLACE_DYADIC<T_GRID>* laplace=(LAPLACE_DYADIC<T_GRID>*)data;
    if(laplace->psi_D(cell->Cell())||laplace->All_Cell_Faces_Neumann(cell))laplace->filled_region_colors(cell->Cell())=-1;
    else laplace->filled_region_colors(cell->Cell())=0;
}
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
Find_Solution_Regions()
{
    FLOOD_FILL flood_fill;filled_region_colors.Resize(grid.number_of_cells);
    // set domain boundary cells and cells with objects to uncolorable (by setting all of the cells to -1, we end up ignoring invalid indices too)
    ARRAYS_COMPUTATIONS::Fill(filled_region_colors,-1);
    MAP_MESH::Map_Cells(grid.uniform_grid,grid.cells,0,this,Mark_Interior_Cells<T_GRID>);
    filled_region_touches_dirichlet.Remove_All();
    // do the fill
    number_of_regions=flood_fill.Flood_Fill(grid,filled_region_colors,psi_N,&filled_region_touches_dirichlet);
}
//#####################################################################
// Function Set_Neumann_Outer_Boundaries
//#####################################################################
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
Set_Neumann_Outer_Boundaries()
{
    for(FACE_ITERATOR iterator(grid,grid.Map_Boundary_Faces());iterator.Valid();iterator.Next())psi_N(iterator.Face_Index())=true;
    pcg.Enforce_Compatibility();
}
//#####################################################################
// Function Set_Dirichlet_Outer_Boundaries
//#####################################################################
template<class T_GRID> static void Set_Boundary_Psi_D(void* data,const typename T_GRID::CELL* cell)
{
    ARRAY<bool>* psi_D=(ARRAY<bool>*)data;
    (*psi_D)(cell->Cell())=true;
}
template<class T_GRID> void LAPLACE_DYADIC<T_GRID>::
Set_Dirichlet_Outer_Boundaries()
{
    MAP_MESH::Map_Ghost_Cells(grid.uniform_grid,grid.cells,1,&psi_D,Set_Boundary_Psi_D<T_GRID>);
}
//#####################################################################
#if COMPILE_WITH_BINTREE_SUPPORT
template class LAPLACE_DYADIC<BINTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_DYADIC<BINTREE_GRID<double> >;
#endif
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class LAPLACE_DYADIC<OCTREE_GRID<float> >;
template class LAPLACE_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_DYADIC<OCTREE_GRID<double> >;
template class LAPLACE_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif

#endif
