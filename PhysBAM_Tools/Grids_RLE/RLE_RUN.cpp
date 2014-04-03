//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_RUN
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_RUN.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_RUN_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_RUN_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Column_Union
//#####################################################################
template<class T_ARRAY_RUN> void RLE_RUN::
Column_Union(const T_ARRAY_RUN& column1,const T_ARRAY_RUN& column2,ARRAY<RLE_RUN>& column_union,const int minimum_long_run_length)
{
    column_union.Remove_All(); 
    const typename T_ARRAY_RUN::ELEMENT *run1=&column1(1),*run2=&column2(1),*run1_end=&column1(column1.m),*run2_end=&column2(column2.m);
    bool current_is_long=false;int new_jmin=run1->jmin;assert(new_jmin==run2->jmin && run1->is_long && run2->is_long);
    while(run1<run1_end && run2<run2_end){
        bool new_is_long=run1->is_long&&run2->is_long;
        if(new_is_long!=current_is_long){current_is_long=new_is_long;
            if(!new_is_long && new_jmin-column_union(column_union.m).jmin<minimum_long_run_length){
                if(column_union.m>1 && !column_union(column_union.m-1).is_long) column_union.m--;
                else column_union(column_union.m).is_long=false;}
            else column_union.Append(RLE_RUN(new_is_long,new_jmin));}
        if((run1+1)->jmin==(run2+1)->jmin){run1++;run2++;new_jmin=run1->jmin;}
        else if((run1+1)->jmin<(run2+1)->jmin){run1++;new_jmin=run1->jmin;}else{run2++;new_jmin=run2->jmin;}}
    column_union.Append(*run1); // append the sentinel
    assert(column_union(1).is_long);
}
//#####################################################################
// Function Dilate_Vertical
//#####################################################################
void RLE_RUN::
Dilate_Vertical(const ARRAY<RLE_RUN>& column,ARRAY<RLE_RUN>& dilated_column,const int extra_cells)
{
    // TODO: write this asserting that long and short runs exactly alternate
    dilated_column.Preallocate(10);int jmin=column(1).jmin,jmax=column(column.m).jmin;
    int last_j=jmin;assert(column(1).is_long);
    for(int r=2;r<column.m;r++)if(!column(r).is_long){
        if(last_j<column(r).jmin-extra_cells){
            dilated_column.Append(RLE_RUN(true,last_j));
            dilated_column.Append(RLE_RUN(false,column(r).jmin-extra_cells));}
        else if(!dilated_column.m) dilated_column.Append(RLE_RUN(false,jmin));
        last_j=column(r+1).jmin+extra_cells;}
    if(last_j<jmax) dilated_column.Append(RLE_RUN(true,last_j));
    dilated_column.Append(RLE_RUN(true,jmax));
    dilated_column.Compact();
}
//#####################################################################
// Function Dilate_Horizontal
//#####################################################################
void RLE_RUN::
Dilate_Horizontal(const ARRAY<RLE_RUN>* columns,ARRAY<RLE_RUN>* new_columns,const int m_start,const int m_end,const int stride,const int extra_cells,const int minimum_long_run_length)
{
    ARRAY<ARRAY<RLE_RUN> > even_unions(extra_cells);
    ARRAY<RLE_RUN> new_pair;new_pair.Preallocate(10);
    Column_Union(columns[stride*m_start],columns[stride*(m_start+1)],new_pair,minimum_long_run_length);
    ARRAYS_COMPUTATIONS::Fill(even_unions,new_pair);
    for(int i=m_start+1;;){
        if(i-extra_cells>=m_start){ // complete left union
            ARRAY<RLE_RUN>& left=new_columns[stride*(i-extra_cells)];
            if(i-2*extra_cells>=m_start){left.Preallocate(10);Column_Union(columns[stride*(i-2*extra_cells)],even_unions(1),left,minimum_long_run_length);left.Compact();}
            else left=even_unions(1);}
        if(++i-extra_cells>m_end) break;
        if(i-extra_cells>=m_start){ // complete right union
            ARRAY<RLE_RUN>& right=new_columns[stride*(i-extra_cells)];
            if(i<=m_end){right.Preallocate(10);Column_Union(even_unions(1),columns[stride*i],right,minimum_long_run_length);right.Compact();}
            else right=even_unions(1);}
        if(++i-extra_cells>m_end) break;
        if(i<=m_end) Column_Union(columns[stride*(i-1)],columns[stride*i],new_pair,minimum_long_run_length);
        else new_pair=columns[stride*m_end];
        for(int u=1;u<extra_cells;u++) Column_Union(new_pair,even_unions(u+1),even_unions(u),minimum_long_run_length);
        ARRAY<RLE_RUN>::Exchange_Arrays(new_pair,even_unions(extra_cells));}
}
//#####################################################################
// Function Split_At_Ground
//#####################################################################
// also converts between ARRAY<RLE_RUN> and ARRAY<T_RUN>
template<class T_RUN> void RLE_RUN::
Split_At_Ground(const ARRAY<RLE_RUN>& column,ARRAY<T_RUN>& new_column,const int ground_j,const int minimum_long_run_length)
{
    assert(column(1).is_long);
    if(ground_j>=column(1).jmin+minimum_long_run_length && ground_j<column(2).jmin){ // split first run at ground
        bool split_run=column(2).jmin-ground_j>=minimum_long_run_length;
        new_column.Resize(column.m+split_run);assert(new_column.m>2);
        new_column(1)=new_column(2)=(T_RUN)column(1);
        for(int r=2;r<=column.m;r++)new_column(r+split_run)=(T_RUN)column(r);
        assert(split_run || !column(2).is_long);new_column(2).jmin=ground_j;}
    else{ // no ground split required
        new_column.Resize(column.m);
        for(int r=1;r<=column.m;r++)new_column(r)=(T_RUN)column(r);}
    if(new_column.m==2){ // hack because face y iterator needs 3 runs minimum in !include_boundary case, TODO: fix
        new_column.Resize(3);
        new_column(3)=new_column(2);
        new_column(2).jmin=(new_column(1).jmin+new_column(2).jmin)/2;}
}
//#####################################################################
template void RLE_RUN::Column_Union(const ARRAY<RLE_RUN_2D>&,const ARRAY<RLE_RUN_2D>&,ARRAY<RLE_RUN>&,const int);
template void RLE_RUN::Column_Union(const ARRAY<RLE_RUN_3D>&,const ARRAY<RLE_RUN_3D>&,ARRAY<RLE_RUN>&,const int);
template void RLE_RUN::Column_Union(const ARRAY<RLE_RUN>&,const ARRAY<RLE_RUN>&,ARRAY<RLE_RUN>&,const int);
template void RLE_RUN::Split_At_Ground(const ARRAY<RLE_RUN>&,ARRAY<RLE_RUN_2D>&,const int,const int);
template void RLE_RUN::Split_At_Ground(const ARRAY<RLE_RUN>&,ARRAY<RLE_RUN_3D>&,const int,const int);
#endif
