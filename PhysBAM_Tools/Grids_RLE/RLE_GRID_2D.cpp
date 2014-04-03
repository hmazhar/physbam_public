//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_2D
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Math_Tools/minabs.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> RLE_GRID_2D<T>::
RLE_GRID_2D()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RLE_GRID_2D<T>::
~RLE_GRID_2D()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Initialize(const RLE_INITIAL_PHI_HELPER<T,2>& initial_phi,const ARRAY<int,VECTOR<int,1> >* ground_j)
{
    Clean_Memory();assert(uniform_grid.counts.x && uniform_grid.counts.y);
    columns.Resize(horizontal_grid.Get_MAC_Grid().Domain_Indices(number_of_ghost_cells+1),false,false);
    int j_start=1-number_of_ghost_cells,j_end=uniform_grid.counts.y+number_of_ghost_cells;
    jmin=j_start-2*minimum_vertical_space;jmax=j_end+2*minimum_vertical_space;
    domain.min_corner.y=uniform_grid.Axis_X(jmin,2);
    ARRAY<RLE_RUN> column;column.Preallocate(10);
    for(int i=columns.domain.min_corner.x;i<=columns.domain.max_corner.x;i++){
        column.Remove_All();
        bool current_is_long=true;int j=j_start;
        column.Append(RLE_RUN_2D(current_is_long,jmin));
        while(j<j_end){
            T phi=initial_phi(uniform_grid.Center(i,j));
            bool new_is_long=phi<-negative_bandwidth||phi>positive_bandwidth;
            if(current_is_long!=new_is_long){
                if(new_is_long) column.Append(RLE_RUN_2D(new_is_long,j));
                else if(column(column.m).jmin>j-minimum_long_run_length) column(column.m).is_long=false;
                else column.Append(RLE_RUN_2D(false,j));
                current_is_long=new_is_long;}
            j+=max(1,(int)(minabs(phi+negative_bandwidth,phi-positive_bandwidth)*uniform_grid.one_over_dX.y));}
        if(!current_is_long) column.Append(RLE_RUN_2D(true,j_end));
        column.Append(RLE_RUN_2D(true,jmax));
        RLE_RUN::Split_At_Ground(column,columns(i),ground_j?(*ground_j)(i):jmin,minimum_long_run_length);}
    if(mpi_grid) mpi_grid->Exchange_Boundary_Columns(number_of_ghost_cells+1,columns);
    Compute_Auxiliary_Information();
}
//#####################################################################
// Function Adjust_Vertical_Space
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Adjust_Vertical_Space()
{
    jmin=1;jmax=uniform_grid.counts.y;
    for(int i=columns.domain.min_corner.x;i<=columns.domain.max_corner.x;i++){
        ARRAY<RLE_RUN_2D>& column=columns(i);assert(column(1).is_long && column(column.m-1).is_long);
        jmin=min(jmin,column(2).jmin-minimum_vertical_space);
        jmax=max(jmax,column(column.m-1).jmin+minimum_vertical_space);}
    if(mpi_grid) mpi_grid->Synchronize_J_Bounds(jmin,jmax);
    for(int i=columns.domain.min_corner.x;i<=columns.domain.max_corner.x;i++){
        ARRAY<RLE_RUN_2D>& column=columns(i);
        column(1).jmin=jmin;column(column.m).jmin=jmax;}
    domain=RANGE<VECTOR<T,2> >(uniform_grid.domain.min_corner.x,uniform_grid.domain.max_corner.x,uniform_grid.Axis_X(jmin,2),uniform_grid.Axis_X(jmax,2));
} 
//#####################################################################
// Function Rebuild
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Rebuild(const RLE_GRID_2D<T>& old_grid,RLE_GRID_2D<T>& new_grid,const ARRAY<bool>& cell_should_be_long,const int extra_cells,const ARRAY<int,VECTOR<int,1> >* ground_j)
{
    int m_start=old_grid.columns.domain.min_corner.x,m_end=old_grid.columns.domain.max_corner.x;int minimum_long_run_length=old_grid.minimum_long_run_length;
    // contract and then dilate each column
    ARRAY<ARRAY<RLE_RUN> ,VECTOR<int,1> > dilated_columns(m_start,m_end,false);
    ARRAY<RLE_RUN> contracted_column;contracted_column.Preallocate(10);
    for(int i=m_start+1;i<=m_end-1;i++){
        // contract
        contracted_column.Remove_All();
        bool current_is_long=false;assert(old_grid.columns(i)(1).is_long);
        for(CELL_ITERATOR cell(old_grid,RANGE<VECTOR<int,1> >(i,i));cell;cell++){ // single column iterator
            bool new_is_long=cell.Long()||cell_should_be_long(cell.Cell());
            if(new_is_long!=current_is_long){current_is_long=new_is_long;contracted_column.Append(RLE_RUN(new_is_long,cell.j));}}
        assert(current_is_long);contracted_column.Append(RLE_RUN(true,old_grid.jmax));
        // dilate
        RLE_RUN::Dilate_Vertical(contracted_column,dilated_columns(i),extra_cells);}

    if(old_grid.mpi_grid) old_grid.mpi_grid->Exchange_Boundary_Columns(extra_cells,dilated_columns);

    // fill in sentinel columns
    dilated_columns(m_start)=dilated_columns(m_start+1);dilated_columns(m_end)=dilated_columns(m_end-1);
    // dilate along x dimension
    if(extra_cells){
        ARRAY<ARRAY<RLE_RUN> ,VECTOR<int,1> > new_dilated_columns(m_start,m_end,false);
        RLE_RUN::Dilate_Horizontal(&dilated_columns(0),&new_dilated_columns(0),m_start,m_end,1,extra_cells,minimum_long_run_length);
        ARRAY<ARRAY<RLE_RUN> ,VECTOR<int,1> >::Exchange_Arrays(dilated_columns,new_dilated_columns);}

    if(old_grid.mpi_grid) old_grid.mpi_grid->Exchange_Boundary_Columns(old_grid.number_of_ghost_cells+1,dilated_columns);

    // construct new_grid
    Transfer_Noncolumn_Data(old_grid,new_grid);
    new_grid.columns.Resize(m_start,m_end,false);
    for(int i=m_start;i<=m_end;i++)
        RLE_RUN::Split_At_Ground(dilated_columns(i),new_grid.columns(i),ground_j?(*ground_j)(i):old_grid.jmin,minimum_long_run_length);
    new_grid.Compute_Auxiliary_Information(); 
}
//#####################################################################
// Function Compute_Auxiliary_Information
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Compute_Auxiliary_Information()
{
    horizontal_grid=uniform_grid.Get_Horizontal_Grid();
    domain=RANGE<VECTOR<T,2> >(uniform_grid.domain.min_corner.x,uniform_grid.domain.max_corner.x,uniform_grid.Axis_X(jmin,2),uniform_grid.Axis_X(jmax,2));
    last_face_numbers.y=number_of_faces=number_of_cells=FACE_Y_ITERATOR::Compute_Indices(*this); // cell numbers are the same as face y numbers
    FACE_X_ITERATOR::Compute_Indices(*this,number_of_faces);last_face_numbers.x=number_of_faces;
    // fill in sentinel indices to make Last_Index_In_Column simpler
    columns(columns.domain.max_corner.x)(1).cell=number_of_cells+1;
    columns(columns.domain.max_corner.x)(1).faces[1]=last_face_numbers.x+1;
    // compute block columns
    Compute_Block_Columns();
}
//#####################################################################
// Function Compute_Block_Columns
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Compute_Block_Columns()
{
    int m_start=columns.domain.min_corner.x,m_end=columns.domain.max_corner.x; 
    // count blocks
    block_run_offsets.Resize(m_start+1,m_end-1,false,false);
    block_run_offsets(m_start+1)=1;
    int index=1;
    for(int i=m_start+1;i<=m_end-2;i++){
        RUN *run1=&columns(i)(1),*run2=&columns(i+1)(1),*run1_end=&columns(i)(columns(i).m);
        while(run1<run1_end){
            int jmax1=(run1+1)->jmin-2,jmax2=(run2+1)->jmin-2;
            int jmin=max(run1->jmin,run2->jmin),jmax=min(jmax1,jmax2);
            if(!run1->is_long && !run2->is_long && jmin<=jmax) index++;
            if(jmax1==jmax) run1++;if(jmax2==jmax) run2++;}
        block_run_offsets(i+1)=index;}
    // compute blocks
    block_runs.Resize(index-1);index=1;
    number_of_blocks=0;
    for(int i=m_start+1;i<=m_end-2;i++){
        RUN *run1=&columns(i)(1),*run2=&columns(i+1)(1),*run1_end=&columns(i)(columns(i).m);
        while(run1<run1_end){
            int jmax1=(run1+1)->jmin-2,jmax2=(run2+1)->jmin-2;
            int jmin=max(run1->jmin,run2->jmin),jmax=min(jmax1,jmax2);
            if(!run1->is_long && !run2->is_long && jmin<=jmax){
                BLOCK_RUN& block=block_runs(index++);block.i=i;block.jmin=jmin;block.jmax=jmax;
                block.block=number_of_blocks+1;number_of_blocks+=jmax-jmin+1;
                int dj1=jmin-run1->jmin,dj2=jmin-run2->jmin;
                block.cells[0]=run1->cell+dj1;block.cells[1]=run2->cell+dj2;
                block.faces_x[0]=run1->faces[0]+dj1;block.faces_x[1]=run1->faces[1]+dj1;block.faces_x[2]=run2->faces[1]+dj2;}
            if(jmax1==jmax) run1++;if(jmax2==jmax) run2++;}}
    assert(index==block_runs.m+1);
}
//#####################################################################
// Function Calculate_Short_Cell_Neighbors
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Calculate_Short_Cell_Neighbors() const
{
    short_cell_neighbors.Resize(number_of_cells);
    for(FACE_X_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell();short_cell_neighbors(c2)(1)=c1;short_cell_neighbors(c1)(2)=c2;}
    for(FACE_Y_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell();short_cell_neighbors(c2)(3)=c1;short_cell_neighbors(c1)(4)=c2;}
}
//#####################################################################
// Function Calculate_Short_Face_Neighbors
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Calculate_Short_Face_Neighbors() const
{
    short_face_neighbors.Resize(number_of_faces);
    for(CELL_ITERATOR cell(*this,number_of_ghost_cells);cell;cell++)if(cell.Short()){
        int fx1=cell.Face_X(0),fx2=cell.Face_X(1);short_face_neighbors(fx2)(1)=fx1;short_face_neighbors(fx1)(2)=fx2; // face x neighbors in x direction
        int fy1=cell.Face_Y(),fy2=fy1+1;short_face_neighbors(fy2)(3)=fy1;short_face_neighbors(fy1)(4)=fy2;} // face y neighbors in y direction
    for(FACE_X_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Short()){
        if(face.cell1.dj || face.cell2.dj){int fx2=face.Face(),fx1=fx2-1;short_face_neighbors(fx2)(3)=fx1;short_face_neighbors(fx1)(4)=fx2;} // face x neighbors in y direction
        if(face.Both_Cells_Short()){
            int fy1=face.cell1.Face_Y(),fy2=face.cell2.Face_Y(),fy3=fy1+1,fy4=fy2+1;
            short_face_neighbors(fy2)(1)=fy1;short_face_neighbors(fy1)(2)=fy2;short_face_neighbors(fy4)(1)=fy3;short_face_neighbors(fy3)(2)=fy4;}} // face y neighbors in x direction
}
//#####################################################################
// Function Calculate_Short_Face_Locations
//#####################################################################
template<class T_GRID> struct Calculate_Short_Face_Locations_Helper{template<class T_FACE> static void Apply(const T_GRID& grid)
{
    for(T_FACE face(grid,grid.number_of_ghost_cells);face;face++)if(face.Short()) grid.short_face_locations(face.Face())=face.X();
}};
template<class T> void RLE_GRID_2D<T>::
Calculate_Short_Face_Locations() const
{
    short_face_locations.Resize(number_of_faces);
    Horizontal_Face_Loop<Calculate_Short_Face_Locations_Helper<RLE_GRID_2D<T> > >(*this);
    for(FACE_Y_ITERATOR face(*this,number_of_ghost_cells);face;face++)short_face_locations(face.Face())=face.X();
}
//#####################################################################
// Function Find_Ghost_Regions
//#####################################################################
template<class T> void RLE_GRID_2D<T>::
Find_Ghost_Regions(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& sentinels,const int number_of_ghost_cells,const bool include_corners) const
{
    RANGE<VECTOR<int,1> > box=horizontal_grid.Get_MAC_Grid().Domain_Indices()+sentinels;
    RANGE<VECTOR<int,1> > band_x(-number_of_ghost_cells,-1);
    regions.Resize(2);
    regions(1)=RANGE<VECTOR<int,1> >(box.min_corner.x,box.min_corner.x)+band_x;
    regions(2)=RANGE<VECTOR<int,1> >(box.max_corner.x,box.max_corner.x)-band_x;
}
//#####################################################################
template class RLE_GRID_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RLE_GRID_2D<double>;
#endif
#endif
