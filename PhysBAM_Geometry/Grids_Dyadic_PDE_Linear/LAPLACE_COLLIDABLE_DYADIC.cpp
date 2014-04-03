//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/FLOOD_FILL_OCTREE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/FLOOD_FILL_QUADTREE.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_NXN.h>
#include <PhysBAM_Geometry/Grids_Dyadic_PDE_Linear/LAPLACE_COLLIDABLE_DYADIC.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_COLLIDABLE_DYADIC<T_GRID>::
Solve(const T time,const bool solution_regions_already_computed)
{
    if(!solution_regions_already_computed) Find_Solution_Regions();

    if(use_second_order_pressure) // propagate all the psi_N values up to the larger faces
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
            if(psi_N(iterator.Face_Index())&&iterator.Deepest_Cell()->Depth_Of_This_Cell()>iterator.Other_Cell()->Depth_Of_This_Cell())
                psi_N(iterator.Other_Face_Index())=true;

    // build all the matrices
    ARRAY<ARRAY<int> > matrix_index_to_cell_index_array(number_of_regions);ARRAY<int> cell_index_to_matrix_index(grid.number_of_cells);
    ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-1,number_of_regions);
    for(int i=1;i<=grid.number_of_cells;i++) filled_region_cell_count(filled_region_colors(i))++;
    ARRAY<SPARSE_MATRIX_NXN<T> > A_array(number_of_regions);ARRAY<VECTOR_ND<T> > b_array(number_of_regions);
    for(int color=1;color<=number_of_regions;color++) if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));
        A_array(color).Set_Size(filled_region_cell_count(color));b_array(color).Resize(filled_region_cell_count(color));}
    ARRAYS_COMPUTATIONS::Fill(filled_region_cell_count,0); // reusing this array in order to make the indirection arrays
    for(int i=1;i<=grid.number_of_cells;i++){
        int color=filled_region_colors(i);
        if(color>0&&(filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            cell_index_to_matrix_index(i)=++filled_region_cell_count(color);
            matrix_index_to_cell_index_array(color)(filled_region_cell_count(color))=i;}}
    Find_A(A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);
    for(int color=1;color<=number_of_regions;color++) if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        pcg.Enforce_Compatibility(!filled_region_touches_dirichlet(color));
        Solve_Subregion(matrix_index_to_cell_index_array(color),A_array(color),b_array(color));}
    if(!solve_neumann_regions) for(int i=1;i<=grid.number_of_cells;i++){
        int filled_region_color=filled_region_colors(i);if(filled_region_color>0 && !filled_region_touches_dirichlet(filled_region_color)) u(i)=0;}
}
//#####################################################################
// Function Solve_Subregion
//#####################################################################
template<class T_GRID> void LAPLACE_COLLIDABLE_DYADIC<T_GRID>::
Solve_Subregion(ARRAY<int>& matrix_index_to_cell_index,SPARSE_MATRIX_NXN<T>& A,VECTOR_ND<T>& b)
{
    int number_of_unknowns=matrix_index_to_cell_index.m;
    A.Negate();b*=(T)-1;
    VECTOR_ND<T> x(number_of_unknowns),q,s,r,k,z;
    for(int i=1;i<=number_of_unknowns;i++) x(i)=u(matrix_index_to_cell_index(i));
    Find_Tolerance(b); // needs to happen after b is completely set up
    if(pcg.show_results) LOG::cout << "solving " << number_of_unknowns << " cells to tolerance " << tolerance << std::endl;
    ARRAY<int> row_lengths(A.n,false);for(int i=1;i<=A.n;i++) row_lengths(i)=A.A(i)->number_of_active_indices;
    SPARSE_MATRIX_FLAT_NXN<T> A_flat;A_flat.Set_Row_Lengths(row_lengths);
    for(int i=1;i<=A.n;i++){
        const SPARSE_VECTOR_ND<T>& vector=*A.A(i);
        for(int index=1;index<=vector.number_of_active_indices;index++) A_flat.Set_Element(i,vector.indices[index],vector.x[index]);}
    pcg.Solve(A_flat,x,b,q,s,r,k,z,tolerance);
    for(int i=1;i<=number_of_unknowns;i++)u(matrix_index_to_cell_index(i))=x(i);
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T_GRID> void LAPLACE_COLLIDABLE_DYADIC<T_GRID>::
Find_A_Off_Diagonal_Helper(FACE_ITERATOR& iterator,ARRAY<SPARSE_MATRIX_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,ARRAY<ARRAY<T> >& row_sum_array)
{
    int face_index=iterator.Face_Index();int axis=iterator.Axis();CELL* cell1=iterator.First_Cell();CELL* cell2=iterator.Second_Cell();
    // if they are different colors, they are disconnected due to a Neumann boundary
    int current_color;
    if(filled_region_colors(cell1->Cell())==-1){
        if(filled_region_colors(cell2->Cell())==-1) return; // both uncolorable
        current_color=filled_region_colors(cell2->Cell());}
    else{
        if(filled_region_colors(cell2->Cell())==-1) current_color=filled_region_colors(cell1->Cell());
        else if(filled_region_colors(cell1->Cell())==filled_region_colors(cell2->Cell()))
            current_color=filled_region_colors(cell1->Cell());
        else return;} // they are different colors (and not dirichlet)
    if(!filled_region_touches_dirichlet(current_color)&&!solve_neumann_regions) return;

    SPARSE_MATRIX_NXN<T>& A=A_array(current_color);VECTOR_ND<T>& b=b_array(current_color);
    ARRAY<T>& row_sum=row_sum_array(current_color);
    if(!psi_N(iterator.Other_Face_Index())){
        const CELL* smaller_cell=iterator.Deepest_Cell();
        const CELL* larger_cell=iterator.Other_Cell();
        T smaller_dx=smaller_cell->DX()[axis+1],larger_dx=larger_cell->DX()[axis+1];
        int smaller_index=cell_index_to_matrix_index(smaller_cell->Cell()),larger_index=cell_index_to_matrix_index(larger_cell->Cell());

        if(use_second_order_pressure&&smaller_cell->Depth_Of_This_Cell()!=larger_cell->Depth_Of_This_Cell()){
            ARRAY<CELL*> all_smaller_cells;
            larger_cell->Get_All_Face_Neighbors(MAP_MESH::face_by_axis[axis][1-iterator.Deepest_Cell_Number()],all_smaller_cells,&grid);
            T distance=0;for(int i=1;i<=all_smaller_cells.m;i++)distance+=all_smaller_cells(i)->Face_Size()*(all_smaller_cells(i)->DX()[axis+1]+larger_dx);
            distance*=(T).5/larger_cell->Face_Size();
            T entry=smaller_cell->Face_Size()/distance;

            for(int i=1;i<=all_smaller_cells.m;i++){
                const CELL* other_smaller_cell=all_smaller_cells(i);
                assert(!psi_D(other_smaller_cell->Cell()));
                if(other_smaller_cell==smaller_cell)continue;
                T entry=(smaller_cell->Face_Size()/larger_cell->Face_Size())*(other_smaller_cell->Face_Size()/distance);
                row_sum(smaller_index)-=entry;
                A.Add_Element(smaller_index,cell_index_to_matrix_index(other_smaller_cell->Cell()),-entry);}

            assert(!psi_D(larger_cell->Cell()));
            row_sum(larger_index)+=entry;A.Add_Element(larger_index,smaller_index,entry);
            row_sum(smaller_index)+=entry;A.Add_Element(smaller_index,larger_index,entry);}
        else{
            T distance=(T).5*(smaller_dx+larger_dx),entry=smaller_cell->Face_Size()/distance;
            if(!psi_D(smaller_cell->Cell()) && !psi_D(larger_cell->Cell())){ // if neither cell is dirichlet
                row_sum(larger_index)+=entry;A.Add_Element(larger_index,smaller_index,entry);
                row_sum(smaller_index)+=entry;A.Add_Element(smaller_index,larger_index,entry);}
            else if(!psi_D(larger_cell->Cell())){
                if(second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface(levelset->phi(larger_cell->Cell()),levelset->phi(smaller_cell->Cell()))){
                    entry/=max(LEVELSET_UTILITIES<T>::Theta(levelset->phi(larger_cell->Cell()),levelset->phi(smaller_cell->Cell())),second_order_cut_cell_threshold); 
                    row_sum(larger_index)+=entry;b(larger_index)-=entry*u_interface(face_index);}
                else{row_sum(larger_index)+=entry;b(larger_index)-=entry*u(smaller_cell->Cell());}}
            else if(!psi_D(smaller_cell->Cell())){
                if(second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface(levelset->phi(smaller_cell->Cell()),levelset->phi(larger_cell->Cell()))){
                    entry/=max(LEVELSET_UTILITIES<T>::Theta(levelset->phi(smaller_cell->Cell()),levelset->phi(larger_cell->Cell())),second_order_cut_cell_threshold); 
                    row_sum(smaller_index)+=entry;b(smaller_index)-=entry*u_interface(face_index);}
                else{row_sum(smaller_index)+=entry;b(smaller_index)-=entry*u(larger_cell->Cell());}}}}
}
namespace PhysBAM{template<class T,class T_GRID> struct FIND_SYSTEM_HELPER{LAPLACE_COLLIDABLE_DYADIC<T_GRID>* laplace;ARRAY<SPARSE_MATRIX_NXN<T> >* A_array;ARRAY<VECTOR_ND<T> >* b_array;
    ARRAY<int>* cell_index_to_matrix_index;ARRAY<ARRAY<T> >* row_sum_array;};}
template<class T,class T_GRID> static void Find_A_Diagonal_And_b_Helper(void* data,const typename T_GRID::CELL* cell)
{
    FIND_SYSTEM_HELPER<T,T_GRID>* helper=(FIND_SYSTEM_HELPER<T,T_GRID>*)data;
    int current_color=helper->laplace->filled_region_colors(cell->Cell());if(current_color<1) return;
    if(!helper->laplace->filled_region_touches_dirichlet(current_color)&&!helper->laplace->solve_neumann_regions) return;
    LAPLACE_COLLIDABLE_DYADIC<T_GRID>* laplace=helper->laplace;SPARSE_MATRIX_NXN<T>& A=(*helper->A_array)(current_color);VECTOR_ND<T>& b=(*helper->b_array)(current_color);
    ARRAY<int>& cell_index_to_matrix_index=*helper->cell_index_to_matrix_index;ARRAY<T>& row_sum=(*helper->row_sum_array)(current_color);
    int cell_index=cell_index_to_matrix_index(cell->Cell());
    b(cell_index)+=cell->Cell_Size()*laplace->f(cell->Cell());A.Set_Element(cell_index,cell_index,-row_sum(cell_index));
}
template<class T_GRID> void LAPLACE_COLLIDABLE_DYADIC<T_GRID>::
Find_A(ARRAY<SPARSE_MATRIX_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index)
{
    ARRAY<ARRAY<T> > row_sum_array(A_array.m);for(int i=1;i<=row_sum_array.m;i++)row_sum_array(i).Resize(filled_region_cell_count(i));
    FIND_SYSTEM_HELPER<T,T_GRID> helper;helper.laplace=this;helper.A_array=&A_array;helper.b_array=&b_array;helper.cell_index_to_matrix_index=&cell_index_to_matrix_index;
    helper.row_sum_array=&row_sum_array;

    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
        Find_A_Off_Diagonal_Helper(iterator,A_array,b_array,cell_index_to_matrix_index,row_sum_array);

    MAP_MESH::Map_Cells(grid.uniform_grid,grid.cells,0,&helper,Find_A_Diagonal_And_b_Helper<T,T_GRID>);
}
//#####################################################################
// Function Find_Solution_Regions
//#####################################################################
template<class T_GRID> static void Mark_Interior_Cells(void* data,const typename T_GRID::CELL* cell)
{
    LAPLACE_COLLIDABLE_DYADIC<T_GRID>* laplace=(LAPLACE_COLLIDABLE_DYADIC<T_GRID>*)data;
    if(laplace->psi_D(cell->Cell())||laplace->All_Cell_Faces_Neumann(cell))laplace->filled_region_colors(cell->Cell())=-1;
    else laplace->filled_region_colors(cell->Cell())=0;
}
template<class T_GRID> void LAPLACE_COLLIDABLE_DYADIC<T_GRID>::
Find_Solution_Regions()
{
    FLOOD_FILL flood_fill;filled_region_colors.Resize(grid.number_of_cells);
    // set domain boundary cells and cells with objects to uncolorable (by setting all of the cells to -1, we end up ignoring invalid indices too)
    filled_region_colors.Fill(-1);
    MAP_MESH::Map_Cells(grid.uniform_grid,grid.cells,0,this,Mark_Interior_Cells<T_GRID>);
    filled_region_touches_dirichlet.Remove_All();
    // do the fill
    number_of_regions=flood_fill.Flood_Fill(grid,filled_region_colors,psi_N,&filled_region_touches_dirichlet);
}
//#####################################################################
// Function Set_Neumann_Outer_Boundaries
//#####################################################################
template<class T_GRID> void LAPLACE_COLLIDABLE_DYADIC<T_GRID>::
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
template<class T_GRID> void LAPLACE_COLLIDABLE_DYADIC<T_GRID>::
Set_Dirichlet_Outer_Boundaries()
{
    MAP_MESH::Map_Ghost_Cells(grid.uniform_grid,grid.cells,1,&psi_D,Set_Boundary_Psi_D<T_GRID>);
}
//#####################################################################

template class LAPLACE_COLLIDABLE_DYADIC<OCTREE_GRID<float> >;
template class LAPLACE_COLLIDABLE_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_COLLIDABLE_DYADIC<OCTREE_GRID<double> >;
template class LAPLACE_COLLIDABLE_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif
