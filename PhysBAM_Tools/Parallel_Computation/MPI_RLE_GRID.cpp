//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#endif
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Function Package_Cell_Data
//#####################################################################
template<class T_GRID> template<class T2> MPI_PACKAGE MPI_RLE_GRID<T_GRID>::
Package_Cell_Data(ARRAY_BASE<T2,ARRAY<T2> >& data,const RANGE<VECTOR<int,1> >& region) const
{
    return MPI_PACKAGE(data,local_grid.Cells_In_Slices(region.min_corner.x,region.max_corner.x));
}
//#####################################################################
// Function Package_Cell_Data
//#####################################################################
template<class T_GRID> template<class T2> MPI_PACKAGE MPI_RLE_GRID<T_GRID>::
Package_Cell_Data(ARRAY_BASE<T2,ARRAY<T2> >& data,const RANGE<VECTOR<int,2> >& region) const
{
    //if(region.min_corner.x==region.max_corner.x)
        return MPI_PACKAGE(data,local_grid.Cells_In_Slice(region.min_corner.x,region.min_corner.y,region.max_corner.y));
#if 0
    int count=region.max_corner.x-region.min_corner.x+1;
    int block_lengths[count],displacements[count];
    for(int index=0;index<count;index++){
        RANGE<VECTOR<int,1> > indices=local_grid.Cells_In_Slice(region.min_corner.x+index,region.min_corner.y,region.max_corner.y);
        block_lengths[index]=indices.max_corner.x-indices.min_corner.x+1;displacements[index]=indices.min_corner.x;}
    return MPI_PACKAGE(data,MPI_UTILITIES::Datatype<T2>().Create_indexed(count,block_lengths,displacements));
#endif
}
//#####################################################################
// Function Package_Face_Data
//#####################################################################
template<class T> template<class TS> MPI_PACKAGE MPI_RLE_GRID<T>::
Package_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const ARRAY<RANGE<VECTOR<int,1> > >& regions) const
{
    // TODO: now that all velocities are in the same array, we could use Create_indexed
    RANGE<VECTOR<int,1> > indices_x=RLE_GRID_2D<T>::FACE_X_ITERATOR::Indices_In_Slices(local_grid,regions(1).min_corner.x,regions(1).max_corner.x),
        indices_y=RLE_GRID_2D<T>::FACE_Y_ITERATOR::Indices_In_Slices(local_grid,regions(2).min_corner.x,regions(2).max_corner.x);
    int block_lengths[2]={indices_x.max_corner.x-indices_x.min_corner.x+1,indices_y.max_corner.x-indices_y.min_corner.x+1};
    MPI::Aint displacements[2]={(MPI::Aint)&data(indices_x.min_corner.x),(MPI::Aint)&data(indices_y.min_corner.x)};
    return MPI_PACKAGE(MPI_UTILITIES::Datatype<T>().Create_hindexed(2,block_lengths,displacements));
}
//#####################################################################
// Function Package_Face_Data
//#####################################################################
template<class T> template<class TS> MPI_PACKAGE MPI_RLE_GRID<T>::
Package_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const ARRAY<RANGE<VECTOR<int,2> > >& regions) const
{   
    // TODO: now that all velocities are in the same array, we could use Create_indexed
    int count=regions(1).max_corner.x-regions(1).min_corner.x+regions(2).max_corner.x-regions(2).min_corner.x+regions(3).max_corner.x-regions(3).min_corner.x+3;
    int block_lengths[count];MPI::Aint displacements[count];
    int index=0;
    for(int i=regions(1).min_corner.x;i<=regions(1).max_corner.x;i++){
        RANGE<VECTOR<int,1> > indices=RLE_GRID_3D<T>::FACE_X_ITERATOR::Indices_In_Slice(local_grid,i,regions(1).min_corner.y,regions(1).max_corner.y);
        block_lengths[index]=indices.max_corner.x-indices.min_corner.x+1;displacements[index]=(MPI::Aint)&data(indices.min_corner.x);index++;}
    for(int i=regions(2).min_corner.x;i<=regions(2).max_corner.x;i++){
        RANGE<VECTOR<int,1> > indices=RLE_GRID_3D<T>::FACE_Y_ITERATOR::Indices_In_Slice(local_grid,i,regions(2).min_corner.y,regions(2).max_corner.y);
        block_lengths[index]=indices.max_corner.x-indices.min_corner.x+1;displacements[index]=(MPI::Aint)&data(indices.min_corner.x);index++;}
    for(int i=regions(3).min_corner.x;i<=regions(3).max_corner.x;i++){
        RANGE<VECTOR<int,1> > indices=RLE_GRID_3D<T>::FACE_Z_ITERATOR::Indices_In_Slice(local_grid,i,regions(3).min_corner.y,regions(3).max_corner.y);
        block_lengths[index]=indices.max_corner.x-indices.min_corner.x+1;displacements[index]=(MPI::Aint)&data(indices.min_corner.x);index++;}
    return MPI_PACKAGE(MPI_UTILITIES::Datatype<T>().Create_hindexed(count,block_lengths,displacements));
}
//#####################################################################
// Function Package_Common_Face_Data
//#####################################################################
template<class T> template<class TS> MPI_PACKAGE MPI_RLE_GRID<T>::
Package_Common_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const int horizontal_axis,const RANGE<VECTOR<int,1> >& region) const
{
    assert(horizontal_axis==1);
    return MPI_PACKAGE(data,RLE_GRID_2D<T>::FACE_X_ITERATOR::Indices_In_Slices(local_grid,region.min_corner.x,region.max_corner.x));
}
//#####################################################################
// Function Package_Common_Face_Data
//#####################################################################
template<class T> template<class TS> MPI_PACKAGE MPI_RLE_GRID<T>::
Package_Common_Face_Data(ARRAY_BASE<TS,ARRAY<TS> >& data,const int horizontal_axis,const RANGE<VECTOR<int,2> >& region) const
{   
    // TODO: now that all velocities are in the same array, we could use Create_indexed
    assert(horizontal_axis==1 || horizontal_axis==2);
    int count=region.max_corner.x-region.min_corner.x+1;
    int block_lengths[count];MPI::Aint displacements[count];
    int index=0;
    for(int i=region.min_corner.x;i<=region.max_corner.x;i++){
        RANGE<VECTOR<int,1> > indices;
        if(horizontal_axis==1) indices=RLE_GRID_3D<T>::FACE_X_ITERATOR::Indices_In_Slice(local_grid,i,region.min_corner.y,region.max_corner.y);
        else indices=RLE_GRID_3D<T>::FACE_Z_ITERATOR::Indices_In_Slice(local_grid,i,region.min_corner.y,region.max_corner.y);
        block_lengths[index]=indices.max_corner.x-indices.min_corner.x+1;displacements[index]=(MPI::Aint)&data(indices.min_corner.x);index++;}
    return MPI_PACKAGE(MPI_UTILITIES::Datatype<T>().Create_hindexed(count,block_lengths,displacements));
}
//#####################################################################
// Function Exchange_Boundary_Columns
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_HORIZONTAL_COLUMN> void MPI_RLE_GRID<T_GRID>::
Exchange_Boundary_Columns(const int bandwidth,T_ARRAYS_HORIZONTAL_COLUMN& columns)
{   
    int tag=Get_Unique_Tag();
    T_BOX_HORIZONTAL_INT sentinels=T_GRID::CELL_ITERATOR::Sentinels();
    // send
    ARRAY<MPI::Request> requests;
    ARRAY<T_BOX_HORIZONTAL_INT> send_regions;Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true);
    ARRAY<ARRAY<char> > buffers(send_regions.m); // must stay alive until after Wait_All
    for(int n=1;n<=send_regions.m;n++)if(all_neighbor_ranks(n)!=MPI::PROC_NULL)
        requests.Append(ISend_Columns(columns,send_regions,n,tag,buffers(n)));
    // probe and receive
    ARRAY<T_BOX_HORIZONTAL_INT> recv_regions;Find_Boundary_Regions(recv_regions,sentinels,false,RANGE<VECTOR<int,1> >(-bandwidth,-1),true);
    for(int message=1;message<=requests.m;message++){
        MPI::Status probe_status;
        comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Columns(columns,recv_regions,tag,probe_status);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function ISend_Columns
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_HORIZONTAL_COLUMN> MPI::Request MPI_RLE_GRID<T_GRID>::
ISend_Columns(const T_ARRAYS_HORIZONTAL_COLUMN& columns,const ARRAY<T_BOX_HORIZONTAL_INT>& regions,const int neighbor,const int tag,ARRAY<char>& buffer) const
{
    int buffer_size=MPI_UTILITIES::Pack_Size<TV_INT>(*comm);
    for(typename T_HORIZONTAL_GRID::CELL_ITERATOR iterator(local_grid.horizontal_grid,regions(neighbor));iterator.Valid();iterator.Next())
        buffer_size+=MPI_UTILITIES::Pack_Size(columns(iterator.Cell_Index()),*comm);
    buffer.Resize(buffer_size);int position=0;
    MPI_UTILITIES::Pack(all_neighbor_directions(neighbor),buffer,position,*comm);
    for(typename T_HORIZONTAL_GRID::CELL_ITERATOR iterator(local_grid.horizontal_grid,regions(neighbor));iterator.Valid();iterator.Next())
        MPI_UTILITIES::Pack(columns(iterator.Cell_Index()),buffer,position,*comm);
    return comm->Isend(&buffer(1),position,MPI::PACKED,all_neighbor_ranks(neighbor),tag);
}
//#####################################################################
// Function Recv_Columns
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_HORIZONTAL_COLUMN> void MPI_RLE_GRID<T_GRID>::
Recv_Columns(T_ARRAYS_HORIZONTAL_COLUMN& columns,const ARRAY<T_BOX_HORIZONTAL_INT>& regions,const int tag,const MPI::Status& probe_status) const
{
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    comm->Recv(&buffer(1),buffer.m,MPI::PACKED,probe_status.Get_source(),tag);
    TV_HORIZONTAL_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,*comm);
    int neighbor=0;all_neighbor_directions.Find(-direction,neighbor);
    for(typename T_HORIZONTAL_GRID::CELL_ITERATOR iterator(local_grid.horizontal_grid,regions(neighbor));iterator.Valid();iterator.Next())
        MPI_UTILITIES::Unpack(columns(iterator.Cell_Index()),buffer,position,*comm);
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID> template<class T_ARRAYS_HORIZONTAL_COLUMN> void MPI_RLE_GRID<T_GRID>::Exchange_Boundary_Columns(const int,T_ARRAYS_HORIZONTAL_COLUMN&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################

#endif

//#####################################################################
#define INSTANTIATION_HELPER_NON_MPI(T_GRID,T,d) \
    template void MPI_RLE_GRID<T_GRID >::Exchange_Boundary_Columns(const int,T_GRID::ARRAYS_HORIZONTAL::REBIND<ARRAY<T_GRID::RUN> >::TYPE&); \
    template void MPI_RLE_GRID<T_GRID >::Exchange_Boundary_Columns(const int,T_GRID::ARRAYS_HORIZONTAL::REBIND<ARRAY<RLE_RUN> >::TYPE&);
#ifdef USE_MPI
#define INSTANTIATION_HELPER_T(T_GRID,T,d) \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Cell_Data(ARRAY_BASE<T,ARRAY<T> >&,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Cell_Data(ARRAY_BASE<T_GRID::VECTOR_T,ARRAY<T_GRID::VECTOR_T > >&,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Cell_Data(ARRAY_BASE<int,ARRAY<int> >&,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Cell_Data(ARRAY_BASE<bool,ARRAY<bool> >&,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Cell_Data(ARRAY_BASE<SYMMETRIC_MATRIX<T,2>,ARRAY<SYMMETRIC_MATRIX<T,2> > >&,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Cell_Data(ARRAY_BASE<SYMMETRIC_MATRIX<T,3>,ARRAY<SYMMETRIC_MATRIX<T,3> > >&,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Face_Data(ARRAY_BASE<T,ARRAY<T> >&,const ARRAY<T_GRID::BOX_HORIZONTAL_INT>&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Face_Data(ARRAY_BASE<bool,ARRAY<bool> >&,const ARRAY<T_GRID::BOX_HORIZONTAL_INT>&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Face_Data(ARRAY_BASE<VECTOR<T,d>,ARRAY<VECTOR<T,d> > >&,const ARRAY<T_GRID::BOX_HORIZONTAL_INT>&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Face_Data(ARRAY_BASE<MATRIX<T,1>,ARRAY<MATRIX<T,1> > >&,const ARRAY<T_GRID::BOX_HORIZONTAL_INT>&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Face_Data(ARRAY_BASE<SYMMETRIC_MATRIX<T,2>,ARRAY<SYMMETRIC_MATRIX<T,2> > >&,const ARRAY<T_GRID::BOX_HORIZONTAL_INT>&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Face_Data(ARRAY_BASE<SYMMETRIC_MATRIX<T,3>,ARRAY<SYMMETRIC_MATRIX<T,3> > >&,const ARRAY<T_GRID::BOX_HORIZONTAL_INT>&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Common_Face_Data(ARRAY_BASE<T,ARRAY<T> >&,const int,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Common_Face_Data(ARRAY_BASE<bool,ARRAY<bool> >&,const int,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Common_Face_Data(ARRAY_BASE<VECTOR<T,d>,ARRAY<VECTOR<T,d> > >&,const int,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Common_Face_Data(ARRAY_BASE<MATRIX<T,1>,ARRAY<MATRIX<T,1> > >&,const int,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Common_Face_Data(ARRAY_BASE<SYMMETRIC_MATRIX<T,2>,ARRAY<SYMMETRIC_MATRIX<T,2> > >&,const int,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    template MPI_PACKAGE MPI_RLE_GRID<T_GRID >::Package_Common_Face_Data(ARRAY_BASE<SYMMETRIC_MATRIX<T,3>,ARRAY<SYMMETRIC_MATRIX<T,3> > >&,const int,const T_GRID::BOX_HORIZONTAL_INT&) const; \
    INSTANTIATION_HELPER_NON_MPI(T_GRID,T,d)
#else
#define INSTANTIATION_HELPER_T(T_GRID,T,d) INSTANTIATION_HELPER_NON_MPI(T_GRID,T,d)
#endif
#define INSTANTIATION_HELPER(T) INSTANTIATION_HELPER_T(RLE_GRID_2D< T >,T,2);INSTANTIATION_HELPER_T(RLE_GRID_3D< T >,T,3);
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
#endif
