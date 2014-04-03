//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void GRID<TV>::
Initialize(const TV_INT& counts_input,const RANGE<TV>& box,const bool MAC_grid)
{
    counts=counts_input;
    domain=box;
    TV_INT effective_counts(counts);
    if(MAC_grid) MAC_offset=(T).5;
    else{MAC_offset=0;effective_counts-=1;}
    if(!TV::dimension) return;
    numbers_of_cells=TV_INT::Componentwise_Max(effective_counts,TV_INT());
    if(effective_counts.Min()>0){
        dX=domain.Edge_Lengths()/TV(effective_counts);
        one_over_dX=(T)1/dX;}
    else{dX=TV();one_over_dX=TV();}
    min_dX=dX.Min();
}
//#####################################################################
// Function Create_Grid_Given_Cell_Size
//#####################################################################
template<class TV> GRID<TV> GRID<TV>::
Create_Grid_Given_Cell_Size(const RANGE<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    assert(cell_size);
    TV_INT cells(TV_INT(ceil((domain.Edge_Lengths())/cell_size))+2*boundary_nodes);
    TV domain_center=domain.Center(),actual_domain_half_size=(T).5*cell_size*TV(cells);
    return GRID<TV>(cells,RANGE<TV>(domain_center-actual_domain_half_size,domain_center+actual_domain_half_size),mac_grid);
}
//#####################################################################
// Function Create_Even_Sized_Grid_Given_Cell_Size
//#####################################################################
template<class TV> GRID<TV> GRID<TV>::
Create_Even_Sized_Grid_Given_Cell_Size(const RANGE<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    assert(cell_size);
    TV edge_lengths(domain.Edge_Lengths());
    TV_INT cells(TV_INT(ceil(edge_lengths/cell_size))+2*boundary_nodes);
    TV cell_sizes;ARRAYS_COMPUTATIONS::Fill(cell_sizes,cell_size);
    for(int i=1;i<=TV::dimension;i++) if(cells(i)&1) cell_sizes(i)=edge_lengths(i)/++cells(i);
    cells+=2*boundary_nodes; // TODO: Why do we do this again?
    TV domain_center=domain.Center(),actual_domain_half_size=(T).5*cell_sizes*TV(cells);
    return GRID<TV>(cells,RANGE<TV>(domain_center-actual_domain_half_size,domain_center+actual_domain_half_size),mac_grid);
}
//#####################################################################
template void GRID<VECTOR<float,0> >::Initialize(VECTOR<int,0> const&,RANGE<VECTOR<float,0> > const&,bool);
template void GRID<VECTOR<float,1> >::Initialize(VECTOR<int,1> const&,RANGE<VECTOR<float,1> > const&,bool);
template void GRID<VECTOR<float,2> >::Initialize(VECTOR<int,2> const&,RANGE<VECTOR<float,2> > const&,bool);
template void GRID<VECTOR<float,3> >::Initialize(VECTOR<int,3> const&,RANGE<VECTOR<float,3> > const&,bool);
template GRID<VECTOR<float,3> > GRID<VECTOR<float,3> >::Create_Grid_Given_Cell_Size(RANGE<VECTOR<float,3> > const&,float,bool,int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void GRID<VECTOR<double,0> >::Initialize(VECTOR<int,0> const&,RANGE<VECTOR<double,0> > const&,bool);
template void GRID<VECTOR<double,1> >::Initialize(VECTOR<int,1> const&,RANGE<VECTOR<double,1> > const&,bool);
template void GRID<VECTOR<double,2> >::Initialize(VECTOR<int,2> const&,RANGE<VECTOR<double,2> > const&,bool);
template void GRID<VECTOR<double,3> >::Initialize(VECTOR<int,3> const&,RANGE<VECTOR<double,3> > const&,bool);
template GRID<VECTOR<double,3> > GRID<VECTOR<double,3> >::Create_Grid_Given_Cell_Size(RANGE<VECTOR<double,3> > const&,double,bool,int);
#endif
