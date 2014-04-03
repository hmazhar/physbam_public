//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_NXN.h>
#include <PhysBAM_Geometry/Grids_Dyadic_PDE_Linear/POISSON_COLLIDABLE_DYADIC.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_beta_And_Add_Jumps_To_b
//#####################################################################
template<class T_GRID> void POISSON_COLLIDABLE_DYADIC<T_GRID>::
Compute_beta_And_Add_Jumps_To_b(const T dt,const T time)
{
    ARRAY<T> phi_ghost(grid.number_of_cells);
    levelset->boundary->Fill_Ghost_Cells_Cell(grid,levelset->phi,phi_ghost,time);
    
    if(use_variable_beta) Find_Variable_beta();else Find_Constant_beta(phi_ghost);
    if(u_jumps) Add_Jump_To_b(phi_ghost);
    if(beta_un_jumps) Add_Derivative_Jump_To_b(phi_ghost);
}
//#####################################################################
// Function Find_Constant_beta
//#####################################################################
// only set up for jump conditons - doesn't work for Dirichlet boundary conditions yet
template<class T_GRID> void POISSON_COLLIDABLE_DYADIC<T_GRID>::
Find_Constant_beta(const ARRAY<T>& phi_ghost)
{      
    // interior
    T half_width=(T).5*number_of_interface_cells*grid.Minimum_Edge_Length();
    if(GFM || smear_beta){
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
            beta_face(iterator.Face_Index())=LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index())),beta_minus,beta_plus,half_width);}
    else{ // smear 1/beta for the delta function method
        T rho_minus=1/beta_minus,rho_plus=1/beta_plus;
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
            beta_face(iterator.Face_Index())=1/LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index())),rho_minus,rho_plus,half_width);}
    if(GFM){ // adjust beta near interface
        T beta_minus_times_beta_plus=beta_minus*beta_plus;
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            if(!LEVELSET_UTILITIES<T>::Interface(phi_ghost(iterator.First_Cell_Index()),phi_ghost(iterator.Second_Cell_Index()))) continue;
            beta_face(iterator.Face_Index())=beta_minus_times_beta_plus/LEVELSET_UTILITIES<T>::Convex_Average(phi_ghost(iterator.First_Cell_Index()),phi_ghost(iterator.Second_Cell_Index()),beta_minus,beta_plus);}}
}
//#####################################################################
// Function Find_Variable_beta
//#####################################################################
// only set up for Dirichlet boundary conditions - doesn't work for jump conditons yet
template<class T_GRID> void POISSON_COLLIDABLE_DYADIC<T_GRID>::
Find_Variable_beta()
{
    if(beta_given_on_faces) return; // beta_face already set
    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
        beta_face(iterator.Face_Index())=(T).5*(variable_beta(iterator.First_Cell_Index())+variable_beta(iterator.Second_Cell_Index()));
}
//#####################################################################
// Function Find_A
//#####################################################################
namespace PhysBAM{template<class T,class T_GRID> struct FIND_SYSTEM_HELPER{POISSON_COLLIDABLE_DYADIC<T_GRID>* poisson;ARRAY<SPARSE_MATRIX_NXN<T> >* A_array;ARRAY<VECTOR_ND<T> >* b_array;
    ARRAY<int>* cell_index_to_matrix_index;ARRAY<ARRAY<T> >* row_sum_array;};}
template<class T_GRID> void POISSON_COLLIDABLE_DYADIC<T_GRID>::
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

    SPARSE_MATRIX_NXN<T>& A=A_array(current_color);VECTOR_ND<T>& b=b_array(current_color);
    ARRAY<T>& row_sum=row_sum_array(current_color);
    if(!psi_N(face_index)){
        const CELL* smaller_cell=iterator.Deepest_Cell();
        const CELL* larger_cell=iterator.Other_Cell();
        T smaller_dx=smaller_cell->DX()[axis+1],larger_dx=larger_cell->DX()[axis+1];
        int smaller_index=cell_index_to_matrix_index(smaller_cell->Cell()),larger_index=cell_index_to_matrix_index(larger_cell->Cell());
        T distance=(T).5*(smaller_dx+larger_dx),entry=beta_face(face_index)*smaller_cell->Face_Size()/distance;
        if(!psi_D(smaller_cell->Cell()) && !psi_D(larger_cell->Cell())){ // if neither cell is dirichlet
            row_sum(larger_index)+=entry;A.Set_Element(larger_index,smaller_index,entry);
            row_sum(smaller_index)+=entry;A.Set_Element(smaller_index,larger_index,entry);}
        else if(!psi_D(larger_cell->Cell())){
            row_sum(larger_index)+=entry;b(larger_index)-=entry*u(smaller_cell->Cell());}
        else if(!psi_D(smaller_cell->Cell())){
            row_sum(smaller_index)+=entry;b(smaller_index)-=entry*u(larger_cell->Cell());}}
}
template<class T,class T_GRID> static void Find_A_Diagonal_And_b_Helper(void* data,const typename T_GRID::CELL* cell)
{
    FIND_SYSTEM_HELPER<T,T_GRID>* helper=(FIND_SYSTEM_HELPER<T,T_GRID>*)data;
    int current_color=helper->poisson->filled_region_colors(cell->Cell());if(current_color<1) return;
    POISSON_COLLIDABLE_DYADIC<T_GRID>* poisson=helper->poisson;SPARSE_MATRIX_NXN<T>& A=(*helper->A_array)(current_color);VECTOR_ND<T>& b=(*helper->b_array)(current_color);
    ARRAY<int>& cell_index_to_matrix_index=*helper->cell_index_to_matrix_index;ARRAY<T>& row_sum=(*helper->row_sum_array)(current_color);
    int cell_index=cell_index_to_matrix_index(cell->Cell());
    b(cell_index)+=cell->Cell_Size()*poisson->f(cell->Cell());A.Set_Element(cell_index,cell_index,-row_sum(cell_index));
}
template<class T_GRID> void POISSON_COLLIDABLE_DYADIC<T_GRID>::
Find_A(ARRAY<SPARSE_MATRIX_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int>& cell_index_to_matrix_index)
{
    ARRAY<ARRAY<T> > row_sum_array(A_array.m);for(int i=1;i<=row_sum_array.m;i++)row_sum_array(i).Resize(filled_region_cell_count(i));
    FIND_SYSTEM_HELPER<T,T_GRID> helper;helper.poisson=this;helper.A_array=&A_array;helper.b_array=&b_array;helper.cell_index_to_matrix_index=&cell_index_to_matrix_index;
    helper.row_sum_array=&row_sum_array;

    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
        Find_A_Off_Diagonal_Helper(iterator,A_array,b_array,cell_index_to_matrix_index,row_sum_array);

    MAP_MESH::Map_Cells(grid.uniform_grid,grid.cells,0,&helper,::Find_A_Diagonal_And_b_Helper<T,T_GRID>);
}
//#####################################################################
// Function Add_Jump_To_b
//#####################################################################
// b is positive still!
template<class T_GRID> void POISSON_COLLIDABLE_DYADIC<T_GRID>::
Add_Jump_To_b(const ARRAY<T>& phi_ghost)
{
    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        int cell1_index=iterator.First_Cell_Index(),cell2_index=iterator.Second_Cell_Index();
        if(!LEVELSET_UTILITIES<T>::Interface(phi_ghost(cell1_index),phi_ghost(cell2_index))||(psi_D(cell1_index)&&psi_D(cell2_index))) continue;
        if(psi_N(iterator.Face_Index())) continue;
        T smaller_dx=iterator.Face_DX(),larger_dx=iterator.Other_Face_DX(),dx=(T).5*(smaller_dx+larger_dx),ratio=sqr(smaller_dx/larger_dx);
        T jump=beta_face(iterator.Face_Index())/dx*LEVELSET_UTILITIES<T>::Average(phi_ghost(cell1_index),u_jump(cell1_index),phi_ghost(cell2_index),u_jump(cell2_index));
        if(iterator.Deepest_Cell_Is_First_Cell()){
            f(cell1_index)-=LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell1_index))*jump/iterator.First_Cell_DX();
            f(cell2_index)-=ratio*LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell2_index))*jump/iterator.Second_Cell_DX();}
        else{
            f(cell1_index)-=ratio*LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell1_index))*jump/iterator.First_Cell_DX();
            f(cell2_index)-=LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell2_index))*jump/iterator.Second_Cell_DX();}}        
}
//#####################################################################
// Function Add_Derivative_Jump_To_b
//#####################################################################
// b is positive still
template<class T_GRID> void POISSON_COLLIDABLE_DYADIC<T_GRID>::
Add_Derivative_Jump_To_b(const ARRAY<T>& phi_ghost)
{
    bool normals_defined=levelset->normals!=0;levelset->Compute_Normals();

    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        int cell1_index=iterator.First_Cell_Index(),cell2_index=iterator.Second_Cell_Index();
        if(!LEVELSET_UTILITIES<T>::Interface(phi_ghost(cell1_index),phi_ghost(cell2_index))||(psi_D(cell1_index)&&psi_D(cell2_index))) continue;
        if(psi_N(iterator.Face_Index())) continue;
        T smaller_dx=iterator.Face_DX(),larger_dx=iterator.Other_Face_DX(),ratio=sqr(smaller_dx/larger_dx);
        T jump=LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell1_index))*beta_face(iterator.Face_Index())*
            LEVELSET_UTILITIES<T>::Average(phi_ghost(cell1_index),beta_un_jump(cell1_index)*(*levelset->normals)(cell1_index)[iterator.Axis()+1],
                                               phi_ghost(cell2_index),beta_un_jump(cell2_index)*(*levelset->normals)(cell2_index)[iterator.Axis()+1]);
        T theta=LEVELSET_UTILITIES<T>::Theta(phi_ghost(cell1_index),phi_ghost(cell2_index));
        if(iterator.Deepest_Cell_Is_First_Cell()){
            f(cell1_index)-=(1-theta)*jump/(iterator.First_Cell_DX()*LEVELSET_UTILITIES<T>::Heaviside(phi_ghost(cell2_index),beta_minus,beta_plus));
            f(cell2_index)-=ratio*theta*jump/(iterator.Second_Cell_DX()*LEVELSET_UTILITIES<T>::Heaviside(phi_ghost(cell1_index),beta_minus,beta_plus));}
        else{
            f(cell1_index)-=ratio*(1-theta)*jump/(iterator.First_Cell_DX()*LEVELSET_UTILITIES<T>::Heaviside(phi_ghost(cell2_index),beta_minus,beta_plus));
            f(cell2_index)-=theta*jump/(iterator.Second_Cell_DX()*LEVELSET_UTILITIES<T>::Heaviside(phi_ghost(cell1_index),beta_minus,beta_plus));}}

    if(!normals_defined){delete levelset->normals;levelset->normals=0;}
}
//#####################################################################
template class POISSON_COLLIDABLE_DYADIC<OCTREE_GRID<float> >;
template class POISSON_COLLIDABLE_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POISSON_COLLIDABLE_DYADIC<OCTREE_GRID<double> >;
template class POISSON_COLLIDABLE_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif
