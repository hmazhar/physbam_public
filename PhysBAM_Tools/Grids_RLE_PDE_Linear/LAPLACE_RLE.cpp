//#####################################################################
// Copyright 2005, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Data_Structures/FLOOD_FILL_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/GRAPH.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE_PDE_Linear/LAPLACE_RLE.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_RLE<T_GRID>::
Solve(const T time,const bool solution_regions_already_computed,const ARRAY<T>* phi_ghost)
{
    if(!solution_regions_already_computed) Find_Solution_Regions();

    // build all the matrices
    ARRAY<ARRAY<int> > matrix_index_to_cell_index_array(number_of_regions);ARRAY<int> cell_index_to_matrix_index(grid.number_of_cells);
    ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-1,number_of_regions);
    for(int c=1;c<=grid.number_of_cells;c++)if(filled_region_colors(c)>0) filled_region_cell_count(filled_region_colors(c))++;
    ARRAY<ARRAY<int> > row_lengths_array(number_of_regions);
    for(int color=1;color<=number_of_regions;color++)if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));
        row_lengths_array(color).Resize(filled_region_cell_count(color));}
    ARRAYS_COMPUTATIONS::Fill(filled_region_cell_count,0); // reusing this array in order to make the indirection arrays
    if(!mpi_grid) for(int c=1;c<=grid.number_of_cells;c++){
        int color=filled_region_colors(c);if(color<1) continue;
        cell_index_to_matrix_index(c)=++filled_region_cell_count(color);
        matrix_index_to_cell_index_array(color)(filled_region_cell_count(color))=c;}
    else laplace_mpi->Find_Matrix_Indices(filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Row_Lengths(row_lengths_array,cell_index_to_matrix_index);
    ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array(number_of_regions);ARRAY<VECTOR_ND<T> > b_array(number_of_regions);
    for(int color=1;color<=number_of_regions;color++)if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        A_array(color).Set_Row_Lengths(row_lengths_array(color));
        b_array(color).Resize(filled_region_cell_count(color));}
    Find_A(A_array,b_array,cell_index_to_matrix_index,phi_ghost);
    for(int color=1;color<=number_of_regions;color++)if(filled_region_cell_count(color)>0 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
        pcg.Enforce_Compatibility(!filled_region_touches_dirichlet(color));
        Solve_Subregion(matrix_index_to_cell_index_array(color),A_array(color),b_array(color),color);}
    if(!solve_neumann_regions) for(int c=1;c<=grid.number_of_cells;c++)if(filled_region_colors(c)==-2) u(c)=0;
}
//#####################################################################
// Function Solve_Subregion
//#####################################################################
template<class T_GRID> void LAPLACE_RLE<T_GRID>::
Solve_Subregion(const ARRAY<int>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int color)
{
    int number_of_unknowns=matrix_index_to_cell_index.m;
    VECTOR_ND<T> x(number_of_unknowns),q,s,r,k,z;
    for(int i=1;i<=number_of_unknowns;i++)x(i)=u(matrix_index_to_cell_index(i));
    Find_Tolerance(b); // needs to happen after b is completely set up
    if(pcg.show_results) LOG::cout<<"solving "<<number_of_unknowns<<" cells to tolerance "<<tolerance<<std::endl;
    if(!mpi_grid) pcg.Solve(A,x,b,q,s,r,k,z,tolerance);
    else laplace_mpi->Solve(A,x,b,q,s,r,k,z,tolerance,color);
    for(int i=1;i<=number_of_unknowns;i++)u(matrix_index_to_cell_index(i))=x(i);
}
//#####################################################################
// Function Find_Row_Lengths
//#####################################################################
template<class T_GRID> void LAPLACE_RLE<T_GRID>::
Find_Row_Lengths(ARRAY<ARRAY<int> >& row_lengths_array,ARRAY<int>& cell_index_to_matrix_index)
{
    T_GRID::template Horizontal_Face_Loop<Horizontal_Row_Length_Contribution>(*this,row_lengths_array,cell_index_to_matrix_index);
    Vertical_Row_Length_Contribution(row_lengths_array,cell_index_to_matrix_index);
    for(int i=1;i<=row_lengths_array.m;i++) row_lengths_array(i)+=1; // diagonal contribution
}
//#####################################################################
// Function Horizontal_Row_Length_Contribution
//#####################################################################
template<class T_GRID> template<class T_FACE> void LAPLACE_RLE<T_GRID>::
Horizontal_Row_Length_Contribution::Apply(const LAPLACE_RLE<T_GRID>& laplace,ARRAY<ARRAY<int> >& row_lengths_array,ARRAY<int>& cell_index_to_matrix_index)
{
    assert(laplace.grid.long_run_cells==2);
    for(T_FACE face(laplace.grid,0);face;face++)if(!laplace.psi_N(face.Face())){
        int c1=face.cell1.Cell(),m1=cell_index_to_matrix_index(c1),m2=cell_index_to_matrix_index(face.cell2.Cell());if(!m1 || !m2) continue;
        ARRAY<int>& row_lengths=row_lengths_array(laplace.filled_region_colors(c1));
        if(face.cell1.Short() && face.cell2.Short()){row_lengths(m1)++;row_lengths(m2)++;}
        else if(face.cell1.Short()){
            int len=face.cell2.length-1,dj=face.cell1.j-face.cell2.j;
            if(dj<len){row_lengths(m1)++;row_lengths(m2)++;}
            if(dj>0){row_lengths(m1)++;row_lengths(m2+1)++;}}
        else if(face.cell2.Short()){
            int len=face.cell1.length-1,dj=face.cell2.j-face.cell1.j;
            if(dj<len){row_lengths(m1)++;row_lengths(m2)++;}
            if(dj>0){row_lengths(m1+1)++;row_lengths(m2)++;}}
        else{
            if(face.Length()>=3){row_lengths(m1)+=2;row_lengths(m1+1)+=2;row_lengths(m2)+=2;row_lengths(m2+1)+=2;} // handle common case
            else{
                if(face.Length()>=2){row_lengths(m1)++;row_lengths(m1+1)++;row_lengths(m2)++;row_lengths(m2+1)++;}
                int jmin1=face.cell1.j,jmin2=face.cell2.j,jmax1=face.cell1.jmax()-1,jmax2=face.cell2.jmax()-1;
                if(jmin1+1<=jmax2-1){row_lengths(m1+1)++;row_lengths(m2)++;}
                if(jmax1-1>=jmin2+1){row_lengths(m1)++;row_lengths(m2+1)++;}}}}
}
//#####################################################################
// Function Vertical_Row_Length_Contribution
//#####################################################################
template<class T_GRID> void LAPLACE_RLE<T_GRID>::
Vertical_Row_Length_Contribution(ARRAY<ARRAY<int> >& row_lengths_array,ARRAY<int>& cell_index_to_matrix_index) const
{
    assert(grid.long_run_cells==2);
    for(FACE_Y_ITERATOR face(grid,0,false);face;face++)if(!psi_N(face.Face())){
        int c2=face.cell2.Cell(),m1=cell_index_to_matrix_index(c2-1),m2=cell_index_to_matrix_index(c2);if(!m1 || !m2) continue;
        ARRAY<int>& row_lengths=row_lengths_array(filled_region_colors(c2));
        row_lengths(m1)++;row_lengths(m2)++;}
    for(CELL_ITERATOR cell(grid,1);cell;cell++)if(cell.Long() && !psi_N(cell.Face_Y()+1)){
        int c=cell.Cell(),m=cell_index_to_matrix_index(c);if(!m) continue;
        ARRAY<int>& row_lengths=row_lengths_array(filled_region_colors(c));
        row_lengths(m)++;row_lengths(m+1)++;}
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T_GRID> void LAPLACE_RLE<T_GRID>::
Find_A(ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,const ARRAY<T>* phi_ghost)
{
    // compute A
    T_GRID::template Horizontal_Face_Loop<Find_Horizontal_Terms>(*this,A_array,b_array,cell_index_to_matrix_index,phi_ghost);
    Find_Vertical_Terms(A_array,b_array,cell_index_to_matrix_index,phi_ghost);
    // compute b
    for(int c=1;c<=grid.number_of_cells;c++){int color=filled_region_colors(c);if(color>=1) b_array(color)(cell_index_to_matrix_index(c))-=f(c);}
}
//#####################################################################
// Function Find_Horizontal_Terms
//#####################################################################
template<class T_GRID> template<class T_FACE> void LAPLACE_RLE<T_GRID>::
Find_Horizontal_Terms::Apply(const LAPLACE_RLE<T_GRID>& laplace,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,
    const ARRAY<T>* phi_ghost)
{
    assert(laplace.grid.long_run_cells==2);
    T area_over_length=laplace.grid.uniform_grid.Face_Size(T_FACE::Axis())/laplace.grid.uniform_grid.dX[T_FACE::Axis()];
    const ARRAY<bool> &psi_N=laplace.psi_N,&psi_D=laplace.psi_D;
    for(T_FACE face(laplace.grid,0);face;face++){int f=face.Face();assert(laplace.grid.long_run_faces_horizontal==1 || face.Short() || psi_N(f)==psi_N(f+1));if(!psi_N(f)){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell(),m1=cell_index_to_matrix_index(c1),m2=cell_index_to_matrix_index(c2);
        int color1=laplace.filled_region_colors(c1),color2=laplace.filled_region_colors(c2),color;
        if(color1>0) color=color1;else if(color2>0) color=color2;else continue;
        SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(color);VECTOR_ND<T>& b=b_array(color);
        if(face.cell1.Short() && psi_D(c2)){
            if(laplace.second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface((*phi_ghost)(c1),(*phi_ghost)(c2))){
                T entry=area_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c1),(*phi_ghost)(c2)),laplace.second_order_cut_cell_threshold);
                A.Add_Element(m1,m1,entry);b(m1)+=entry*laplace.u_interface(f);}
            else{A.Add_Element(m1,m1,area_over_length);b(m1)+=area_over_length*laplace.u(c2);}}
        else if(face.cell2.Short() && psi_D(c1)){
            if(laplace.second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface((*phi_ghost)(c2),(*phi_ghost)(c1))){
                T entry=area_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c2),(*phi_ghost)(c1)),laplace.second_order_cut_cell_threshold);
                A.Add_Element(m2,m2,entry);b(m2)+=entry*laplace.u_interface(f);}
            else{A.Add_Element(m2,m2,area_over_length);b(m2)+=area_over_length*laplace.u(c1);}}
        else if(face.cell1.Short() && face.cell2.Short()){A.Add_Element(m1,m1,area_over_length);A.Add_Element(m2,m2,area_over_length);A.Set_Symmetric_Elements(m1,m2,-area_over_length);}
        else if(face.cell1.Short()){
            assert(!psi_D(c2));
            if(laplace.psi_N(face.cell2.Face_Y()+1)) PHYSBAM_FATAL_ERROR("Invalid long internal psi_N");
            int len=face.cell2.length-1,dj=face.cell1.j-face.cell2.j;
            T a12_hi=area_over_length*dj/len,a12_lo=area_over_length-a12_hi,a22_lo_hi=a12_hi*(len-dj)/len;
            A.Add_Element(m2,m2,a12_lo-a22_lo_hi);A.Add_Element(m2+1,m2+1,a12_hi-a22_lo_hi);A.Add_Symmetric_Elements(m2,m2+1,a22_lo_hi);
            if(psi_D(c1)){
                b(m2)+=a12_lo*laplace.u(c1);
                b(m2+1)+=a12_hi*laplace.u(c1);}
            else{
                A.Add_Element(m1,m1,area_over_length);
                if(dj<len) A.Set_Symmetric_Elements(m1,m2,-a12_lo);
                if(dj>0) A.Set_Symmetric_Elements(m1,m2+1,-a12_hi);}}
        else if(face.cell2.Short()){
            assert(!psi_D(c1));
            if(laplace.psi_N(face.cell1.Face_Y()+1)) PHYSBAM_FATAL_ERROR("Invalid long internal psi_N");
            int len=face.cell1.length-1,dj=face.cell2.j-face.cell1.j;
            T a12_hi=area_over_length*dj/len,a12_lo=area_over_length-a12_hi,a11_lo_hi=a12_hi*(len-dj)/len;
            A.Add_Element(m1,m1,a12_lo-a11_lo_hi);A.Add_Element(m1+1,m1+1,a12_hi-a11_lo_hi);A.Add_Symmetric_Elements(m1,m1+1,a11_lo_hi);
            if(psi_D(c2)){
                b(m1)+=a12_lo*laplace.u(c2);
                b(m1+1)+=a12_hi*laplace.u(c2);}
            else{
                A.Add_Element(m2,m2,area_over_length);
                if(dj<len) A.Set_Symmetric_Elements(m1,m2,-a12_lo);
                if(dj>0) A.Set_Symmetric_Elements(m1+1,m2,-a12_hi);}}
        else{ // TODO: refer the unfortunate reader to some documentation
            if(laplace.psi_N(face.cell1.Face_Y()+1) || laplace.psi_N(face.cell2.Face_Y()+1)) PHYSBAM_FATAL_ERROR("Invalid long internal psi_N");
            long long len_f=face.Length()-1,len_c1=face.cell1.length-1,len_c2=face.cell2.length-1,dj1=face.cell1.j-face.j(),dj2=face.cell2.j-face.j();
            long long a12_sum_sum_scaled=6*(1+len_f)*len_c1*len_c2;
            long long a12_sum_hi_scaled=3*(1+len_f)*len_c1*(len_f-2*dj2);
            long long a12_hi_sum_scaled=3*(1+len_f)*len_c2*(len_f-2*dj1);
            long long a12_hi_hi_scaled=(1+len_f)*(len_f*(2*len_f+1-3*(dj1+dj2))+6*dj1*dj2);
            long long a12_lo_hi_scaled=a12_sum_hi_scaled-a12_hi_hi_scaled,a12_hi_lo_scaled=a12_hi_sum_scaled-a12_hi_hi_scaled;
            long long a12_lo_lo_scaled=a12_sum_sum_scaled-a12_sum_hi_scaled-a12_hi_lo_scaled;
            assert(a12_lo_lo_scaled==(1+len_f)*(len_f*(2*len_f+1-3*(len_c2+dj2))-3*(len_c1+dj1)*(len_f-2*(len_c2+dj2))));
            T scale12=area_over_length/(T)(6*len_c1*len_c2),scale11=area_over_length/(T)(6*len_c1*len_c1),scale22=area_over_length/(T)(6*len_c2*len_c2);
            T a12_lo_lo=scale12*(T)a12_lo_lo_scaled,a12_lo_hi=scale12*(T)a12_lo_hi_scaled,a12_hi_lo=scale12*(T)a12_hi_lo_scaled,a12_hi_hi=scale12*(T)a12_hi_hi_scaled;
            T a11_lo_hi=scale11*(T)((1+len_f)*(len_f*(3*len_c1+6*dj1-2*len_f-1)-6*dj1*(len_c1+dj1)));
            T a22_lo_hi=scale22*(T)((1+len_f)*(len_f*(3*len_c2+6*dj2-2*len_f-1)-6*dj2*(len_c2+dj2)));
            if(!psi_D(c1)){A.Add_Element(m1,m1,a12_lo_lo+a12_lo_hi-a11_lo_hi);A.Add_Element(m1+1,m1+1,a12_hi_hi+a12_hi_lo-a11_lo_hi);A.Add_Symmetric_Elements(m1,m1+1,a11_lo_hi);}
            if(!psi_D(c2)){A.Add_Element(m2,m2,a12_lo_lo+a12_hi_lo-a22_lo_hi);A.Add_Element(m2+1,m2+1,a12_hi_hi+a12_lo_hi-a22_lo_hi);A.Add_Symmetric_Elements(m2,m2+1,a22_lo_hi);}
            if(psi_D(c2)){
                b(m1)+=a12_lo_lo*laplace.u(c2)+a12_lo_hi*laplace.u(c2+1);
                b(m1+1)+=a12_hi_lo*laplace.u(c2)+a12_hi_hi*laplace.u(c2+1);}
            else if(psi_D(c1)){
                b(m2)+=a12_lo_lo*laplace.u(c1)+a12_hi_lo*laplace.u(c1+1);
                b(m2+1)+=a12_lo_hi*laplace.u(c1)+a12_hi_hi*laplace.u(c1+1);}
            else{
                if(a12_lo_lo) A.Set_Symmetric_Elements(m1,m2,-a12_lo_lo);
                if(a12_lo_hi) A.Set_Symmetric_Elements(m1,m2+1,-a12_lo_hi);
                if(a12_hi_lo) A.Set_Symmetric_Elements(m1+1,m2,-a12_hi_lo);
                if(a12_hi_hi) A.Set_Symmetric_Elements(m1+1,m2+1,-a12_hi_hi);}}}}
}
//#####################################################################
// Function Find_Vertical_Terms
//#####################################################################
template<class T_GRID> void LAPLACE_RLE<T_GRID>::
Find_Vertical_Terms(ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,ARRAY<int>& cell_index_to_matrix_index,const ARRAY<T>* phi_ghost) const
{
    assert(grid.long_run_cells==2);
    T area_over_length=grid.uniform_grid.Face_Size(2)/grid.uniform_grid.dX.y;
    for(FACE_Y_ITERATOR face(grid,0,false);face;face++){int f=face.Face();if(!psi_N(f)){
        int c2=face.cell2.Cell(),c1=c2-1,m1=cell_index_to_matrix_index(c1),m2=cell_index_to_matrix_index(c2);
        int color1=filled_region_colors(c1),color2=filled_region_colors(c2),color;
        if(color1>0) color=color1;else if(color2>0) color=color2;else continue;
        SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(color);VECTOR_ND<T>& b=b_array(color);
        if(psi_D(c2)){
            if(second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface((*phi_ghost)(c1),(*phi_ghost)(c2))){
                T entry=area_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c1),(*phi_ghost)(c2)),second_order_cut_cell_threshold);
                A.Add_Element(m1,m1,entry);b(m1)+=entry*u_interface(f);}
            else{A.Add_Element(m1,m1,area_over_length);b(m1)+=area_over_length*u(c2);}}
        else if(psi_D(c1)){
            if(second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface((*phi_ghost)(c2),(*phi_ghost)(c1))){
                T entry=area_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c2),(*phi_ghost)(c1)),second_order_cut_cell_threshold);
                A.Add_Element(m2,m2,entry);b(m2)+=entry*u_interface(f);}
            else{A.Add_Element(m2,m2,area_over_length);b(m2)+=area_over_length*u(c1);}}
        else{A.Add_Element(m1,m1,area_over_length);A.Add_Element(m2,m2,area_over_length);A.Set_Symmetric_Elements(m1,m2,-area_over_length);}}}
    for(CELL_ITERATOR cell(grid,1);cell;cell++)if(cell.Long() && !psi_N(cell.Face_Y()+1)){
        int c=cell.Cell(),m=cell_index_to_matrix_index(c),color=filled_region_colors(c);if(color<1) continue;
        SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(color);
        T entry=area_over_length/(cell.length-1);A.Add_Element(m,m,entry);A.Add_Element(m+1,m+1,entry);A.Add_Symmetric_Elements(m,m+1,-entry);}
}
//#####################################################################
// Function Find_Solution_Regions
//#####################################################################
template<class T_GRID> void LAPLACE_RLE<T_GRID>::
Find_Solution_Regions()
{
    assert(grid.long_run_cells==2);
    // setup cell graph
    GRAPH graph(grid.number_of_cells);
    T_GRID::template Horizontal_Face_Loop<Add_Horizontal_Edges>(*this,graph);
    for(FACE_Y_ITERATOR face(grid,0,false);face;face++){int f=face.Face(),c2=face.cell2.Cell(); // don't add vertical edges in ghost layer
        if(!psi_N(f)) graph.Add_Undirected_Edge(c2-1,c2);
        if(face.cell2.Long() && !psi_N(f+1)) graph.Add_Undirected_Edge(c2,c2+1);}
    // set domain boundary cells and cells with objects to uncolorable (by setting all of the cells to -1, we end up ignoring invalid indices too)
    filled_region_colors.Resize(grid.number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(filled_region_colors,-1);
    for(int c=1;c<=grid.number_of_cells;c++)if(graph.edges(c).m && !psi_D(c)) filled_region_colors(c)=0;
    filled_region_touches_dirichlet.Remove_All();
    // do the fill
    FLOOD_FILL_GRAPH flood_fill;
    number_of_regions=flood_fill.Flood_Fill(graph,filled_region_colors,&filled_region_touches_dirichlet);
    // correct flood fill for distributed grids
    if(mpi_grid) laplace_mpi->Synchronize_Solution_Regions();
    // set unsolved neumann regions to color -2
    if(!solve_neumann_regions) for(int c=1;c<=grid.number_of_cells;c++){
        int color=filled_region_colors(c);if(color>0 && !filled_region_touches_dirichlet(color)) filled_region_colors(c)=-2;}
}
//#####################################################################
// Function Add_Horizontal_Edges
//#####################################################################
template<class T_GRID> template<class T_FACE> void LAPLACE_RLE<T_GRID>::
Add_Horizontal_Edges::Apply(const LAPLACE_RLE<T_GRID>& laplace,GRAPH& graph)
{
    for(T_FACE face(laplace.grid,0);face;face++)if(!laplace.psi_N(face.Face())){int c1=face.cell1.Cell(),c2=face.cell2.Cell();
        for(int i1=0;i1<=(int)face.cell1.Long();i1++)for(int i2=0;i2<=(int)face.cell2.Long();i2++)graph.Add_Undirected_Edge(c1+i1,c2+i2);}
}
//#####################################################################
template class LAPLACE_RLE<RLE_GRID_2D<float> >;
template class LAPLACE_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_RLE<RLE_GRID_2D<double> >;
template class LAPLACE_RLE<RLE_GRID_3D<double> >;
#endif
#endif
