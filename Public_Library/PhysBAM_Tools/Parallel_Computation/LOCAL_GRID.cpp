//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parallel_Computation/LOCAL_GRID.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LOCAL_GRID<T_GRID>::
LOCAL_GRID(const T_GRID& global_grid_input)
    :mpi_grid(grid,3,true),global_grid(global_grid_input)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LOCAL_GRID<T_GRID>::
LOCAL_GRID(const T_GRID& global_grid_input,const T_GRID& local_grid_input)
    :grid(local_grid_input),mpi_grid(grid,3,true),global_grid(global_grid_input)
{
    Initialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LOCAL_GRID<T_GRID>::
~LOCAL_GRID()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T_GRID> void LOCAL_GRID<T_GRID>::
Initialize()
{// find offset
    TV X=grid.Domain().Minimum_Corner();
    TV_INT local_offset=grid.Closest_Node(X),global_offset=global_grid.Closest_Node(X);
    offset=global_offset-local_offset;
    PHYSBAM_ASSERT((grid.Node(local_offset)-global_grid.Node(global_offset)).Magnitude()<(T).01*grid.Minimum_Edge_Length(),"mismatch between global and local grids");
    // pretend there are neighbors on all sides for use in Find_Boundary_Regions
    mpi_grid.side_neighbor_ranks.Resize(T_GRID::number_of_neighbors_per_cell);
    ARRAYS_COMPUTATIONS::Fill(mpi_grid.side_neighbor_ranks,1);
    mpi_grid.all_neighbor_ranks.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAYS_COMPUTATIONS::Fill(mpi_grid.all_neighbor_ranks,1);
    ARRAY<RANGE<TV_INT> > regions;
    mpi_grid.Find_Boundary_Regions(regions,RANGE<TV_INT>::Zero_Box(),true,RANGE<VECTOR<int,1> >(-5,-5),true);
    RANGE<TV_INT> global_region=global_grid.Domain_Indices();
    neighbor_overlaps.Resize(regions.m);
    for(int n=1;n<=regions.m;n++) if(regions(n).Lazy_Intersection(global_region-offset)) neighbor_overlaps(n)=true;
}
//#####################################################################
// Function Put
//#####################################################################
template<class T_GRID> template<class T_ARRAYS> void LOCAL_GRID<T_GRID>::
Put(const T_ARRAYS& local_data,const RANGE<TV_INT>& region,T_ARRAYS& global_data) const
{
    CELL_ITERATOR local(grid,region),global(grid,region+offset);
    for(;local.Valid();local.Next(),global.Next()) if(global_data.Valid_Index(global.Cell_Index()) && local_data.Valid_Index(local.Cell_Index()))
        global_data(global.Cell_Index())=local_data(local.Cell_Index());
}
//#####################################################################
// Function Put
//#####################################################################
template<class T_GRID> template<class T_ARRAYS> void LOCAL_GRID<T_GRID>::
Put(const T_ARRAYS& local_data,T_ARRAYS& global_data,const RANGE<TV_INT>& sentinels) const
{
    Put(local_data,Interior_Region(sentinels),global_data);
    ARRAY<RANGE<TV_INT> > regions;
    mpi_grid.Find_Boundary_Regions(regions,sentinels,true,RANGE<VECTOR<int,1> >(-3,-1),true);
    for(int n=1;n<=regions.m;n++) if(!neighbor_overlaps(n)) Put(local_data,regions(n),global_data);
}
//#####################################################################
// Function Get
//#####################################################################
template<class T_GRID> template<class T_ARRAYS> void LOCAL_GRID<T_GRID>::
Get(const T_ARRAYS& global_data,T_ARRAYS& local_data) const
{
    RANGE<TV_INT> region=local_data.Domain_Indices();
    CELL_ITERATOR local(grid,region),global(grid,region+offset);
    for(;local.Valid();local.Next(),global.Next()) local_data(local.Cell_Index())=global_data(global.Cell_Index());
}
//#####################################################################
// Function Put_Faces
//#####################################################################
template<class T_GRID> template<class T_FACE_ARRAYS> void LOCAL_GRID<T_GRID>::
Put_Faces(const T_FACE_ARRAYS& local_data,T_FACE_ARRAYS& global_data) const
{
    for(int axis=1;axis<=T_GRID::dimension;axis++) Put(local_data.Component(axis),global_data.Component(axis),mpi_grid.Face_Sentinels(axis));
}
//#####################################################################
// Function Maximum_Error
//#####################################################################
template<class T_GRID> template<class T_ARRAYS> typename T_GRID::VECTOR_T::SCALAR LOCAL_GRID<T_GRID>::
Maximum_Error(const T_ARRAYS& local_data,const T_ARRAYS& global_data,const int bandwidth,TV_INT& index,const RANGE<TV_INT>& sentinels) const
{
    RANGE<TV_INT> region=Interior_Region(sentinels).Thickened(bandwidth);
    T max_error=0;CELL_ITERATOR local(grid,region),global(grid,region+offset);
    for(;local.Valid();local.Next(),global.Next()) if(global_data.Valid_Index(global.Cell_Index()) && local_data.Valid_Index(local.Cell_Index())){
            T error=sqrt(ARRAYS_COMPUTATIONS::Magnitude_Squared(global_data(global.Cell_Index())-local_data(local.Cell_Index())));
        if(max_error<error){max_error=error;index=local.Cell_Index();}}
    return max_error;
}
//#####################################################################
// Function Maximum_Error
//#####################################################################
template<class T_GRID> template<class T_FACE_ARRAYS> typename T_GRID::VECTOR_T::SCALAR LOCAL_GRID<T_GRID>::
Maximum_Error(const std::string& prefix,const T_FACE_ARRAYS& local_data,const T_FACE_ARRAYS& global_data,const int bandwidth,const T threshold)
{
    T max_error=(T)0;
    const char *axis_names[3]={"x","y","z"};
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        TV_INT index;
        max_error=Maximum_Error(local_data.Component(axis),global_data.Component(axis),bandwidth,index,mpi_grid.Face_Sentinels(axis));
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(max_error>threshold){LOG::cout<<prefix<<", face "<<axis_names[axis-1]<<" = "<<max_error<<" ("<<index<<")"<<std::endl;}
#endif
        }
    return max_error;
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class LOCAL_GRID<GRID<VECTOR<T,d> > >; \
    template T LOCAL_GRID<GRID<VECTOR<T,d> > >::Maximum_Error<ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<GRID<VECTOR<T,d> > >::Maximum_Error<ARRAY<VECTOR<T,d>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d>,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<GRID<VECTOR<T,d> > >::Maximum_Error<ARRAY<T,VECTOR<int,d> > >(ARRAY<T,VECTOR<int,d> > const&,ARRAY<T,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<GRID<VECTOR<T,d> > >::Maximum_Error<ARRAY<bool,VECTOR<int,d> > >(ARRAY<bool,VECTOR<int,d> > const&,ARRAY<bool,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<GRID<VECTOR<T,d> > >::Maximum_Error<ARRAY<int,VECTOR<int,d> > >(ARRAY<int,VECTOR<int,d> > const&,ARRAY<int,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<GRID<VECTOR<T,d> > >::Maximum_Error<ARRAY<bool,FACE_INDEX<d> > >(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,ARRAY<bool,FACE_INDEX<d> > const&,ARRAY<bool,FACE_INDEX<d> > const&,int,T); \
    template T LOCAL_GRID<GRID<VECTOR<T,d> > >::Maximum_Error<ARRAY<T,FACE_INDEX<d> > >(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,ARRAY<T,FACE_INDEX<d> > const&,ARRAY<T,FACE_INDEX<d> > const&,int,T); \
    template void LOCAL_GRID<GRID<VECTOR<T,d> > >::Put<ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<GRID<VECTOR<T,d> > >::Put<ARRAY<VECTOR<T,d>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d>,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<GRID<VECTOR<T,d> > >::Put<ARRAY<T,VECTOR<int,d> > >(ARRAY<T,VECTOR<int,d> > const&,ARRAY<T,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<GRID<VECTOR<T,d> > >::Put<ARRAY<bool,VECTOR<int,d> > >(ARRAY<bool,VECTOR<int,d> > const&,ARRAY<bool,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<GRID<VECTOR<T,d> > >::Put<ARRAY<int,VECTOR<int,d> > >(ARRAY<int,VECTOR<int,d> > const&,ARRAY<int,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<GRID<VECTOR<T,d> > >::Put_Faces<ARRAY<bool,FACE_INDEX<d> > >(ARRAY<bool,FACE_INDEX<d> > const&,ARRAY<bool,FACE_INDEX<d> >&) const; \
    template void LOCAL_GRID<GRID<VECTOR<T,d> > >::Put_Faces<ARRAY<T,FACE_INDEX<d> > >(ARRAY<T,FACE_INDEX<d> > const&,ARRAY<T,FACE_INDEX<d> >&) const;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
