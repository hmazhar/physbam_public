//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Math_Tools/minabs.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RLE_GRID_3D<T>::
RLE_GRID_3D()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RLE_GRID_3D<T>::
~RLE_GRID_3D()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Initialize(const RLE_INITIAL_PHI_HELPER<T,3>& initial_phi,const ARRAY<int,VECTOR<int,2> >* ground_j,const bool verbose)
{
    Clean_Memory();assert(uniform_grid.counts.x && uniform_grid.counts.y && uniform_grid.counts.z);
    columns.Resize(horizontal_grid.Get_MAC_Grid().Domain_Indices(number_of_ghost_cells+1),false,false);
    int j_start=1-number_of_ghost_cells,j_end=uniform_grid.counts.y+number_of_ghost_cells;
    jmin=j_start-2*minimum_vertical_space;jmax=j_end+2*minimum_vertical_space;
    domain.min_corner.y=uniform_grid.Axis_X(jmin,2);
    ARRAY<RLE_RUN> column;column.Preallocate(10);
    PROGRESS_INDICATOR progress(columns.array.Size());
    for(int i=columns.domain.min_corner.x;i<=columns.domain.max_corner.x;i++)for(int ij=columns.domain.min_corner.y;ij<=columns.domain.max_corner.y;ij++){
        column.Remove_All();
        bool current_is_long=true;int j=j_start;
        column.Append(RLE_RUN_3D(current_is_long,jmin));
        while(j<j_end){
            T phi=initial_phi(uniform_grid.Center(i,j,ij));
            bool new_is_long=phi<-negative_bandwidth||phi>positive_bandwidth;
            if(current_is_long!=new_is_long){
                if(new_is_long) column.Append(RLE_RUN_3D(new_is_long,j));
                else if(column(column.m).jmin>j-minimum_long_run_length) column(column.m).is_long=false;
                else column.Append(RLE_RUN_3D(false,j));
                current_is_long=new_is_long;}
            j+=max(1,(int)(minabs(phi+negative_bandwidth,phi-positive_bandwidth)*uniform_grid.one_over_dX.y));}
        if(!current_is_long) column.Append(RLE_RUN_3D(true,j_end));
        column.Append(RLE_RUN_3D(true,jmax));
        RLE_RUN::Split_At_Ground(column,columns(i,ij),ground_j?(*ground_j)(i,ij):jmin,minimum_long_run_length);
        if(verbose) progress.Progress();}
    if(mpi_grid) mpi_grid->Exchange_Boundary_Columns(number_of_ghost_cells+1,columns);
    Compute_Auxiliary_Information();
}
//#####################################################################
// Function Adjust_Vertical_Space
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Adjust_Vertical_Space()
{
    Adjust_Vertical_Space(1,uniform_grid.counts.y);
}
//#####################################################################
// Function Adjust_Vertical_Space
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Adjust_Vertical_Space(const int max_jmin,const int min_jmax)
{
    jmin=max_jmin;jmax=min_jmax;
    for(int i=columns.domain.min_corner.x;i<=columns.domain.max_corner.x;i++)for(int ij=columns.domain.min_corner.y;ij<=columns.domain.max_corner.y;ij++){
        ARRAY<RLE_RUN_3D>& column=columns(i,ij);assert(column(1).is_long && column(column.m-1).is_long);
        jmin=min(jmin,column(2).jmin-minimum_vertical_space);
        jmax=max(jmax,column(column.m-1).jmin+minimum_vertical_space);}
    if(mpi_grid) mpi_grid->Synchronize_J_Bounds(jmin,jmax);
    for(int i=columns.domain.min_corner.x;i<=columns.domain.max_corner.x;i++)for(int ij=columns.domain.min_corner.y;ij<=columns.domain.max_corner.y;ij++){
        ARRAY<RLE_RUN_3D>& column=columns(i,ij);
        column(1).jmin=jmin;column(column.m).jmin=jmax;}
    domain=RANGE<TV>(uniform_grid.domain.min_corner.x,uniform_grid.domain.max_corner.x,uniform_grid.Axis_X(jmin,2),uniform_grid.Axis_X(jmax,2),uniform_grid.domain.min_corner.z,uniform_grid.domain.max_corner.z);
}
//#####################################################################
// Function Rebuild
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Rebuild(const RLE_GRID_3D<T>& old_grid,RLE_GRID_3D<T>& new_grid,const ARRAY<bool>& cell_should_be_long,const int extra_cells,const ARRAY<int,VECTOR<int,2> >* ground_j)
{
    int m_start=old_grid.columns.domain.min_corner.x,m_end=old_grid.columns.domain.max_corner.x,n_start=old_grid.columns.domain.min_corner.y,n_end=old_grid.columns.domain.max_corner.y;int minimum_long_run_length=old_grid.minimum_long_run_length;
    // contract and then dilate each column
    ARRAY<ARRAY<RLE_RUN> ,VECTOR<int,2> > dilated_columns(m_start,m_end,n_start,n_end,false);
    ARRAY<RLE_RUN> contracted_column;contracted_column.Preallocate(10);
    for(int i=m_start+1;i<=m_end-1;i++)for(int ij=n_start+1;ij<=n_end-1;ij++){
        // contract
        contracted_column.Remove_All();
        bool current_is_long=false;assert(old_grid.columns(i,ij)(1).is_long);
        for(CELL_ITERATOR cell(old_grid,RANGE<VECTOR<int,2> >(i,i,ij,ij));cell;cell++){ // single column iterator
            bool new_is_long=cell.Long()||cell_should_be_long(cell.Cell());
            if(new_is_long!=current_is_long){current_is_long=new_is_long;contracted_column.Append(RLE_RUN(new_is_long,cell.j));}}
        assert(current_is_long);contracted_column.Append(RLE_RUN(true,old_grid.jmax));
        // dilate
        RLE_RUN::Dilate_Vertical(contracted_column,dilated_columns(i,ij),extra_cells);}

    if(old_grid.mpi_grid) old_grid.mpi_grid->Exchange_Boundary_Columns(extra_cells,dilated_columns);

    // fill in sentinel columns
    for(int ij=n_start+1;ij<=n_end-1;ij++){dilated_columns(m_start,ij)=dilated_columns(m_start+1,ij);dilated_columns(m_end,ij)=dilated_columns(m_end-1,ij);}
    for(int i=m_start;i<=m_end;i++){dilated_columns(i,n_start)=dilated_columns(i,n_start+1);dilated_columns(i,n_end)=dilated_columns(i,n_end-1);}
    // dilate along x and z dimensions
    if(extra_cells){
        ARRAY<ARRAY<RLE_RUN> ,VECTOR<int,2> > new_dilated_columns(m_start,m_end,n_start,n_end,false);
        for(int ij=n_start;ij<=n_end;ij++)RLE_RUN::Dilate_Horizontal(&dilated_columns(0,ij),&new_dilated_columns(0,ij),m_start,m_end,dilated_columns.counts.y,extra_cells,minimum_long_run_length);
        for(int i=m_start;i<=m_end;i++)RLE_RUN::Dilate_Horizontal(&new_dilated_columns(i,0),&dilated_columns(i,0),n_start,n_end,1,extra_cells,minimum_long_run_length);}

    if(old_grid.mpi_grid) old_grid.mpi_grid->Exchange_Boundary_Columns(old_grid.number_of_ghost_cells+1,dilated_columns);

    // construct new_grid
    Transfer_Noncolumn_Data(old_grid,new_grid);
    new_grid.columns.Resize(m_start,m_end,n_start,n_end,false);
    for(int i=m_start;i<=m_end;i++)for(int ij=n_start;ij<=n_end;ij++)
        RLE_RUN::Split_At_Ground(dilated_columns(i,ij),new_grid.columns(i,ij),ground_j?(*ground_j)(i,ij):old_grid.jmin,minimum_long_run_length);
    new_grid.Compute_Auxiliary_Information();
}
//#####################################################################
// Function Compute_Auxiliary_Information
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Compute_Auxiliary_Information()
{
    horizontal_grid=uniform_grid.Get_Horizontal_Grid();
    domain=RANGE<TV>(uniform_grid.domain.min_corner.x,uniform_grid.domain.max_corner.x,uniform_grid.Axis_X(jmin,2),uniform_grid.Axis_X(jmax,2),uniform_grid.domain.min_corner.z,uniform_grid.domain.max_corner.z);
    last_face_numbers.y=number_of_faces=number_of_cells=FACE_Y_ITERATOR::Compute_Indices(*this); // cell numbers are the same as face y numbers
    FACE_X_ITERATOR::Compute_Indices(*this,number_of_faces);last_face_numbers.x=number_of_faces;
    FACE_Z_ITERATOR::Compute_Indices(*this,number_of_faces);last_face_numbers.z=number_of_faces;
    // fill in sentinel indices to make Last_Index_In_Column simpler
    columns(columns.domain.max_corner.x,columns.domain.min_corner.y+1)(1).cell=number_of_cells+1;
    columns(columns.domain.max_corner.x,columns.domain.min_corner.y+1)(1).faces[4]=last_face_numbers.z+1;
    for(int i=columns.domain.min_corner.x+1;i<columns.domain.max_corner.x;i++){
        RLE_RUN_3D& run=columns(i,columns.domain.max_corner.y)(1),next_run=columns(i+1,columns.domain.min_corner.y+1)(1);
        run.cell=next_run.cell;run.faces[5]=next_run.faces[4];run.faces[0]=next_run.faces[0];}
    columns(columns.domain.max_corner.x,columns.domain.max_corner.y)(1).faces[0]=last_face_numbers.x+1;
    // compute block columns
    Compute_Block_Columns();
}
//#####################################################################
// Function Compute_Block_Columns
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Compute_Block_Columns()
{
    int m_start=columns.domain.min_corner.x,m_end=columns.domain.max_corner.x,n_start=columns.domain.min_corner.y,n_end=columns.domain.max_corner.y;
    // count blocks
    block_run_offsets.Resize(m_start+1,m_end-2,n_start+1,n_end-1,false,false);
    int index=1;
    for(int i=m_start+1;i<=m_end-2;i++){
        block_run_offsets(i,n_start+1)=index;
        for(int ij=n_start+1;ij<=n_end-2;ij++){
            RUN *run1=&columns(i,ij)(1),*run2=&columns(i+1,ij)(1),*run3=&columns(i,ij+1)(1),*run4=&columns(i+1,ij+1)(1),*run1_end=&columns(i,ij)(columns(i,ij).m);
            while(run1<run1_end){
                int jmax1=(run1+1)->jmin-2,jmax2=(run2+1)->jmin-2,jmax3=(run3+1)->jmin-2,jmax4=(run4+1)->jmin-2;
                int jmin=max(run1->jmin,run2->jmin,run3->jmin,run4->jmin),jmax=min(jmax1,jmax2,jmax3,jmax4);
                if(!run1->is_long && !run2->is_long && !run3->is_long && !run4->is_long && jmin<=jmax) index++;
                if(jmax1==jmax) run1++;if(jmax2==jmax) run2++;if(jmax3==jmax) run3++;if(jmax4==jmax) run4++;}
            block_run_offsets(i,ij+1)=index;}}
    // compute blocks
    block_runs.Resize(index-1);index=1;
    number_of_blocks=0;
    for(int i=m_start+1;i<=m_end-2;i++)
        for(int ij=n_start+1;ij<=n_end-2;ij++){
            RUN *run1=&columns(i,ij)(1),*run2=&columns(i+1,ij)(1),*run3=&columns(i,ij+1)(1),*run4=&columns(i+1,ij+1)(1),*run1_end=&columns(i,ij)(columns(i,ij).m);
            while(run1<run1_end){
                int jmax1=(run1+1)->jmin-2,jmax2=(run2+1)->jmin-2,jmax3=(run3+1)->jmin-2,jmax4=(run4+1)->jmin-2;
                int jmin=max(run1->jmin,run2->jmin,run3->jmin,run4->jmin),jmax=min(jmax1,jmax2,jmax3,jmax4);
                if(!run1->is_long && !run2->is_long && !run3->is_long && !run4->is_long && jmin<=jmax){
                    BLOCK_RUN& block=block_runs(index++);block.ij=ij;block.jmin=jmin;block.jmax=jmax;
                    block.block=number_of_blocks+1;number_of_blocks+=jmax-jmin+1;
                    int dj1=jmin-run1->jmin,dj2=jmin-run2->jmin,dj3=jmin-run3->jmin,dj4=jmin-run4->jmin;
                    block.cells[0]=run1->cell+dj1;block.cells[1]=run2->cell+dj2;block.cells[2]=run3->cell+dj3;block.cells[3]=run4->cell+dj4;
                    block.faces_x[0]=run1->faces[0]+dj1;block.faces_x[1]=run1->faces[1]+dj1;block.faces_x[2]=run2->faces[1]+dj2;
                    block.faces_x[3]=run3->faces[0]+dj3;block.faces_x[4]=run3->faces[1]+dj3;block.faces_x[5]=run4->faces[1]+dj4;
                    block.faces_z[0]=run1->faces[4]+dj1;block.faces_z[1]=run2->faces[4]+dj2;block.faces_z[2]=run1->faces[5]+dj1;
                    block.faces_z[3]=run2->faces[5]+dj2;block.faces_z[4]=run3->faces[5]+dj3;block.faces_z[5]=run4->faces[5]+dj4;}
                if(jmax1==jmax) run1++;if(jmax2==jmax) run2++;if(jmax3==jmax) run3++;if(jmax4==jmax) run4++;}}
    assert(index==block_runs.m+1);
}
//#####################################################################
// Function Calculate_Short_Cell_Neighbors
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Calculate_Short_Cell_Neighbors() const
{
    short_cell_neighbors.Resize(number_of_cells);
    for(FACE_X_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell();short_cell_neighbors(c2)(1)=c1;short_cell_neighbors(c1)(2)=c2;}
    for(FACE_Z_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell();short_cell_neighbors(c2)(5)=c1;short_cell_neighbors(c1)(6)=c2;}
    for(FACE_Y_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell();short_cell_neighbors(c2)(3)=c1;short_cell_neighbors(c1)(4)=c2;}
}
//#####################################################################
// Function Calculate_Short_Face_Neighbors
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Calculate_Short_Face_Neighbors() const
{
    short_face_neighbors.Resize(number_of_faces);
    for(CELL_ITERATOR cell(*this,number_of_ghost_cells);cell;cell++)if(cell.Short()){
        int fx1=cell.Face_X(0),fx2=cell.Face_X(1);short_face_neighbors(fx2)(1)=fx1;short_face_neighbors(fx1)(2)=fx2; // face x neighbors in x direction
        int fz1=cell.Face_Z(0),fz2=cell.Face_Z(1);short_face_neighbors(fz2)(5)=fz1;short_face_neighbors(fz1)(6)=fz2; // face z neighbors in z direction
        int fy1=cell.Face_Y(),fy2=fy1+1;short_face_neighbors(fy2)(3)=fy1;short_face_neighbors(fy1)(4)=fy2;} // face y neighbors in y direction
    for(FACE_X_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Short()){
        if(face.cell1.dj || face.cell2.dj){int fx2=face.Face(),fx1=fx2-1;short_face_neighbors(fx2)(3)=fx1;short_face_neighbors(fx1)(4)=fx2;} // face x neighbors in y direction
        if(face.Both_Cells_Short()){
            int fy1=face.cell1.Face_Y(),fy2=face.cell2.Face_Y(),fy3=fy1+1,fy4=fy2+1;
            short_face_neighbors(fy2)(1)=fy1;short_face_neighbors(fy1)(2)=fy2;short_face_neighbors(fy4)(1)=fy3;short_face_neighbors(fy3)(2)=fy4; // face y neighbors in x direction
            int fz1=face.cell1.Face_Z(0),fz2=face.cell2.Face_Z(0),fz3=face.cell1.Face_Z(1),fz4=face.cell2.Face_Z(1);
            short_face_neighbors(fz2)(1)=fz1;short_face_neighbors(fz1)(2)=fz2;short_face_neighbors(fz4)(1)=fz3;short_face_neighbors(fz3)(2)=fz4;}} // face z neighbors in x direction
    for(FACE_Z_ITERATOR face(*this,number_of_ghost_cells,false);face;face++)if(face.Short()){
        if(face.cell1.dj || face.cell2.dj){int fz2=face.Face(),fz1=fz2-1;short_face_neighbors(fz2)(3)=fz1;short_face_neighbors(fz1)(4)=fz2;} // face z neighbors in y direction
        if(face.Both_Cells_Short()){
            int fy1=face.cell1.Face_Y(),fy2=face.cell2.Face_Y(),fy3=fy1+1,fy4=fy2+1;
            short_face_neighbors(fy2)(5)=fy1;short_face_neighbors(fy1)(6)=fy2;short_face_neighbors(fy4)(5)=fy3;short_face_neighbors(fy3)(6)=fy4; // face y neighbors in z direction
            int fx1=face.cell1.Face_X(0),fx2=face.cell2.Face_X(0),fx3=face.cell1.Face_X(1),fx4=face.cell2.Face_X(1);
            short_face_neighbors(fx2)(5)=fx1;short_face_neighbors(fx1)(6)=fx2;short_face_neighbors(fx4)(5)=fx3;short_face_neighbors(fx3)(6)=fx4;}} // face x neighbors in z direction
}
//#####################################################################
// Function Calculate_Short_Face_Locations
//#####################################################################
template<class T_GRID> struct Calculate_Short_Face_Locations_Helper{template<class T_FACE> static void Apply(const T_GRID& grid)
{
    for(T_FACE face(grid,grid.number_of_ghost_cells);face;face++)if(face.Short())
        grid.short_face_locations(face.Face())=face.X();
}};
template<class T> void RLE_GRID_3D<T>::
Calculate_Short_Face_Locations() const
{
    short_face_locations.Resize(number_of_faces);
    Horizontal_Face_Loop<Calculate_Short_Face_Locations_Helper<RLE_GRID_3D<T> > >(*this);
    for(FACE_Y_ITERATOR face(*this,number_of_ghost_cells);face;face++)short_face_locations(face.Face())=face.X();

}
//#####################################################################
// Function Find_Ghost_Regions
//#####################################################################
template<class T> void RLE_GRID_3D<T>::
Find_Ghost_Regions(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& sentinels,const int number_of_ghost_cells,const bool include_corners,const bool use_separate_regions_for_corners) const
{
    RANGE<VECTOR<int,2> > box=horizontal_grid.Get_MAC_Grid().Domain_Indices()+sentinels; 
    RANGE<VECTOR<int,1> > band(-number_of_ghost_cells,-1);
    RANGE<VECTOR<int,2> > band_x(band.min_corner.x,band.max_corner.x,0,0),band_z(0,0,band.min_corner.x,band.max_corner.x);
    regions.Resize(include_corners&&use_separate_regions_for_corners?8:4);
    // sides
    regions(1)=RANGE<VECTOR<int,2> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.max_corner.y)+band_x;
    regions(2)=RANGE<VECTOR<int,2> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y)-band_x;
    regions(3)=RANGE<VECTOR<int,2> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y)+band_z;
    regions(4)=RANGE<VECTOR<int,2> >(box.min_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y)-band_z;
    // corners
    if(include_corners){
        if(!use_separate_regions_for_corners){
            RANGE<VECTOR<int,2> > extend(-number_of_ghost_cells,number_of_ghost_cells,0,0);
            regions(3)+=extend;
            regions(4)+=extend;}
        else{
            regions(5)=VECTOR<int,2>(box.min_corner.x,box.min_corner.y)+band_x+band_z;
            regions(6)=VECTOR<int,2>(box.min_corner.x,box.max_corner.y)+band_x-band_z;
            regions(7)=VECTOR<int,2>(box.max_corner.x,box.min_corner.y)-band_x+band_z;
            regions(8)=VECTOR<int,2>(box.max_corner.x,box.max_corner.y)-band_x-band_z;}}
}
//#####################################################################
template class RLE_GRID_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RLE_GRID_3D<double>;
#endif
}
#endif
