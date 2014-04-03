//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/Hash.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> THREADED_UNIFORM_GRID<T_GRID>::
THREADED_UNIFORM_GRID(ARRAY<THREAD_PACKAGE>& buffers_input,const int tid_input,const int number_of_threads,T_GRID& local_grid_input,const int number_of_ghost_cells_input,
    const bool skip_initialization,const TV_INT& processes_per_dimension,const TV_BOOL& periodic_input)
    :MPI_GRID<T_GRID>(local_grid_input,number_of_ghost_cells_input,true,processes_per_dimension,periodic_input,0),tid(tid_input),buffers(buffers_input)
{
    number_of_processes=number_of_threads;rank=tid-1;
    if(skip_initialization) return;

    if(tid==1) LOG::SCOPE scope("THREADED INITIALIZE","Initializing Threading");
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(tid==1) LOG::cout<<"number of processes = "<<number_of_processes<<std::endl;
#endif

#ifdef USE_PTHREADS
    if(tid==1){
        lock=new pthread_mutex_t;
        pthread_mutex_init(lock,NULL);
        THREAD_PACKAGE pack(sizeof(pthread_mutex_t*));
        *(pthread_mutex_t**)(&pack.buffer(1))=lock;
        buffers.Append(pack);}
    else{
        while(buffers.m==0) sleep(1);
        lock=*(pthread_mutex_t**)(&buffers(1).buffer(1));}
    if(tid==1){
        barr=new pthread_barrier_t;
        pthread_barrier_init(barr,NULL,number_of_processes);
        THREAD_PACKAGE pack(sizeof(pthread_barrier_t*));
        *(pthread_barrier_t**)(&pack.buffer(1))=barr;
        buffers.Append(pack);}
    else{
        while(buffers.m<=1) sleep(1);
        barr=*(pthread_barrier_t**)(&buffers(2).buffer(1));}
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif

    // extract global grid and divide among processes
    global_grid=local_grid.Get_Regular_Grid(); // Split_Grid and Initialize_Communicator currently assume grid is non-MAC
    Split_Grid(processes_per_dimension);

    // setup topology
    process_ranks.Resize(process_grid.Domain_Indices(1));ARRAYS_COMPUTATIONS::Fill(process_ranks.array,-1);
    TV_INT extents=process_grid.Domain_Indices().Maximum_Corner();
    // sort axes in decreasing order of how much we have to communicate along them
    ARRAY<int> axes(T_GRID::dimension);ARRAY<T> axis_lengths(T_GRID::dimension);
    for(int axis=1;axis<=T_GRID::dimension;axis++){axes(axis)=axis;axis_lengths(axis)=(T)global_grid.Domain_Indices().Maximum_Corner()[axis]/extents[axis];}
    Sort(axes,Indirect_Comparison(axis_lengths));
    // lay out process ranks on grid
    Fill_Process_Ranks(process_grid,process_ranks,axes);
    // fill in ghost process_ranks for periodic domains
    if(periodic!=TV_BOOL()) for(NODE_ITERATOR iterator(process_grid,1,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        TV_INT node=iterator.Node_Index(),wrapped_node=node;
        for(int axis=1;axis<=T_GRID::dimension;axis++) if(periodic[axis]) wrapped_node[axis]=(node[axis]+process_grid.Counts()[axis]-1)%process_grid.Counts()[axis]+1;
        process_ranks(node)=process_ranks(wrapped_node);}
    all_coordinates.Resize(number_of_processes);
    for(NODE_ITERATOR iterator(process_grid);iterator.Valid();iterator.Next())
        all_coordinates(process_ranks(iterator.Node_Index())+1)=iterator.Node_Index();
    coordinates=all_coordinates(tid);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout<<"process_ranks = \n"<<process_ranks<<std::endl;
    LOG::cout<<"coordinates = "<<coordinates<<std::endl;
#endif
    
    side_neighbor_ranks.Resize(T_GRID::number_of_neighbors_per_node);
    side_neighbor_directions.Resize(T_GRID::number_of_neighbors_per_node);
    all_neighbor_ranks.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    all_neighbor_directions.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++){
        side_neighbor_ranks(n)=process_ranks(T_GRID::Node_Neighbor(coordinates,n));
        side_neighbor_directions(n)=T_GRID::Node_Neighbor(TV_INT(),n);}
    for(int n=1;n<=T_GRID::number_of_one_ring_neighbors_per_cell;n++){
        all_neighbor_ranks(n)=process_ranks(T_GRID::One_Ring_Neighbor(coordinates,n));
        all_neighbor_directions(n)=T_GRID::One_Ring_Neighbor(TV_INT(),n);}

    // restrict this process to correct piece
    local_grid=Restrict_Grid(coordinates);
    // initialize offset
    TV_INT start_index;TV_INT end_index;
    for(int axis=1;axis<=T_GRID::dimension;axis++)
        start_index[axis]=boundaries(axis)(coordinates[axis]);
    local_to_global_offset=start_index-TV_INT::All_Ones_Vector();
    // initialize global column index boundaries
    global_column_index_boundaries.Resize(number_of_processes);
    int offset=0;
    for(int proc=1;proc<=all_coordinates.m;proc++){
        TV_INT proc_coordinates=all_coordinates(proc);
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            start_index[axis]=boundaries(axis)(proc_coordinates[axis]);
            end_index[axis]=boundaries(axis)(proc_coordinates[axis]+1);}
        int owned_pressure_values=(end_index-start_index).Product();
        global_column_index_boundaries(proc)=VECTOR<int,2>(offset+1,offset+owned_pressure_values);
        offset+=owned_pressure_values;}
    // initialize the column indices of the ring of values that will be in the domain as well as one cell around it
    local_cell_index_to_global_column_index_map.Resize(local_grid.Domain_Indices(1));
    // first go through the interior and set what those n values will be
    offset=global_column_index_boundaries(tid).x;
    for(CELL_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next()){
        local_cell_index_to_global_column_index_map(iterator.Cell_Index())=offset;
        offset++;}
    // now go through each of the boundaries and do those
    for(int axis=1;axis<=T_GRID::dimension;axis++)
        for(int axis_side=1;axis_side<=2;axis_side++){
            int side=2*axis-(2-axis_side);
            int neighbor_rank=side_neighbor_ranks(side);
            if(neighbor_rank>=0){
                TV_INT axis_vector=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT::Axis_Vector(axis);
                int start_column_index=global_column_index_boundaries(neighbor_rank+1).x;
                // Make that neighbors local_grid
                TV_INT proc_coordinates=all_coordinates(neighbor_rank+1);
                for(int temp_axis=1;temp_axis<=T_GRID::dimension;temp_axis++){
                    start_index[temp_axis]=boundaries(temp_axis)(proc_coordinates[temp_axis]);
                    end_index[temp_axis]=boundaries(temp_axis)(proc_coordinates[temp_axis]+1);}
                CELL_ITERATOR my_iterator(local_grid,0,T_GRID::BOUNDARY_INTERIOR_REGION,side);
                int neighbor_side=axis_side==1?2*axis:2*axis-1;
                T_GRID neighbor_grid=T_GRID(end_index-start_index+TV_INT::All_Ones_Vector(),RANGE<TV>(global_grid.X(start_index),global_grid.X(end_index))).Get_MAC_Grid();
                for(CELL_ITERATOR neighbor_iterator(neighbor_grid,0,T_GRID::BOUNDARY_INTERIOR_REGION,neighbor_side);neighbor_iterator.Valid();neighbor_iterator.Next(),my_iterator.Next()){
                    TV_INT my_cell_index=my_iterator.Cell_Index();
                    local_cell_index_to_global_column_index_map(my_cell_index+axis_vector)=neighbor_iterator.Flat_Index()+start_column_index-1;}}}
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,1> >& process_grid,ARRAY<int,VECTOR<int,1> >& process_ranks,ARRAY<int>& axes)
{
    VECTOR<int,1> extents=process_grid.Domain_Indices().Maximum_Corner();
    for(int i=1;i<=extents.x;i++)process_ranks(i)=i-1;
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,2> >& process_grid,ARRAY<int,VECTOR<int,2> >& process_ranks,ARRAY<int>& axes)
{
    int next_rank=0;
    VECTOR<int,2> extents=process_grid.Domain_Indices().Maximum_Corner(),half_extents=extents/2;
    for(int i=1;i<=half_extents[axes(1)];i++)for(int j=1;j<=half_extents[axes(2)];j++)
        for(int ii=0;ii<2;ii++)for(int jj=0;jj<2;jj++){
            VECTOR<int,2> permuted_index(2*i+ii-1,2*j+jj-1),index;
            for(int a=1;a<=2;a++)index[axes(a)]=permuted_index[a];
            process_ranks(index)=next_rank++;}
    for(int i=1;i<=extents.x;i++)for(int j=1;j<=extents.y;j++)if(process_ranks(i,j)==-1) process_ranks(i,j)=next_rank++;
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,3> >& process_grid,ARRAY<int,VECTOR<int,3> >& process_ranks,ARRAY<int>& axes)
{
    int next_rank=0;
    VECTOR<int,3> extents=process_grid.Domain_Indices().Maximum_Corner(),half_extents=extents/2;
    for(int i=1;i<=half_extents[axes(1)];i++)for(int j=1;j<=half_extents[axes(2)];j++)for(int ij=1;ij<=half_extents[axes(3)];ij++)
        for(int ii=0;ii<2;ii++)for(int jj=0;jj<2;jj++)for(int ijij=0;ijij<2;ijij++){
            VECTOR<int,3> permuted_index(2*i+ii-1,2*j+jj-1,2*ij+ijij-1),index;
            for(int a=1;a<=3;a++)index[axes(a)]=permuted_index[a];
            process_ranks(index)=next_rank++;}
    for(int i=1;i<=extents.x;i++)for(int j=1;j<=extents.y;j++)for(int ij=1;ij<=extents.z;ij++)if(process_ranks(i,j,ij)==-1) process_ranks(i,j,ij)=next_rank++;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T_GRID> void THREADED_UNIFORM_GRID<T_GRID>::
Initialize(VECTOR<VECTOR<bool,2>,T_GRID::dimension>& domain_walls)
{   
    // fix walls
    for(int i=1;i<=T_GRID::number_of_neighbors_per_node;i++)
        if(side_neighbor_ranks(i)!=-1) domain_walls((i+1)/2)(i&1?1:2)=false;

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout<<"rank = "<<rank<<std::endl;
    for(int axis=1;axis<=T_GRID::dimension;axis++)
        LOG::cout<<"boundaries "<<axis<<" = "<<boundaries(axis)(coordinates[axis])<<" to "<<boundaries(axis)(coordinates[axis]+1)<<std::endl;
    LOG::cout<<"topology = "<<process_grid.Domain_Indices()<<std::endl;
    LOG::cout<<"process ranks = \n"<<process_ranks;
#endif
}
//#####################################################################
// Function Synchronize_Dt
//#####################################################################
template<class T_GRID> void THREADED_UNIFORM_GRID<T_GRID>::
Synchronize_Dt(T& dt) const
{
#ifdef USE_PTHREADS
    THREAD_PACKAGE pack(sizeof(T));*(T*)(&pack.buffer(1))=dt;
    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers.m;i++) dt=min(*(T*)(&buffers(i).buffer(1)),dt);
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function All_Reduce
//#####################################################################
template<class T_GRID> void THREADED_UNIFORM_GRID<T_GRID>::
All_Reduce(bool& flag) const
{
#ifdef USE_PTHREADS
    THREAD_PACKAGE pack(sizeof(bool));*(bool*)(&pack.buffer(1))=flag;
    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers.m;i++) flag|=*(bool*)(&buffers(i).buffer(1));
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Exchange_Boundary_Cell_Data
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Exchange_Boundary_Cell_Data(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const int bandwidth,const bool include_corners) const
{
#ifdef USE_PTHREADS
    RANGE<TV_INT> sentinels=RANGE<TV_INT>::Zero_Box();
    const ARRAY<int>& neighbor_ranks=include_corners?all_neighbor_ranks:side_neighbor_ranks;
    // send
    ARRAY<RANGE<TV_INT> > send_regions;Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),include_corners,true,local_grid);
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=-1){
        THREAD_PACKAGE pack=Package_Cell_Data(data,send_regions(n));pack.recv_tid=neighbor_ranks(n);
        pthread_mutex_lock(lock);
        buffers.Append(pack);
        pthread_mutex_unlock(lock);}
    // receive
    ARRAY<RANGE<TV_INT> > recv_regions;Find_Boundary_Regions(recv_regions,sentinels,false,RANGE<VECTOR<int,1> >(-bandwidth,-1),include_corners,true,local_grid);
    pthread_barrier_wait(barr);
    for(int n=1;n<=recv_regions.m;n++)if(neighbor_ranks(n)!=-1){int index=0;
        for(int i=1;i<=buffers.m;i++) if(buffers(i).send_tid==neighbor_ranks(n) && buffers(i).recv_tid==rank) index=i;
        PHYSBAM_ASSERT(index);int position=0;
        for(CELL_ITERATOR iterator(local_grid,recv_regions(n));iterator.Valid();iterator.Next()) data.Unpack(buffers(index).buffer,position,iterator.Cell_Index());}
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Exchange_Boundary_Face_Data
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Exchange_Boundary_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data,const int bandwidth) const
{   
#ifdef USE_PTHREADS
    RANGE<VECTOR<int,1> > boundary_band(0,bandwidth-1),ghost_band(-bandwidth,-1);
    // send
    ARRAY<ARRAY<RANGE<TV_INT> > > send_regions(T_GRID::dimension);
    for(int axis=1;axis<=T_GRID::dimension;axis++)Find_Boundary_Regions(send_regions(axis),Face_Sentinels(axis),true,boundary_band,true);
    for(int n=1;n<=send_regions(1).m;n++)if(all_neighbor_ranks(n)!=-1){
        ARRAY<RANGE<TV_INT> > send_regions_n(T_GRID::dimension);for(int axis=1;axis<=T_GRID::dimension;axis++)send_regions_n(axis)=send_regions(axis)(n);
        int size=0;for(int axis=1;axis<=T_GRID::dimension;axis++) for(FACE_ITERATOR iterator(local_grid,send_regions_n(axis),axis);iterator.Valid();iterator.Next()) size+=data.data(axis).Pack_Size();
        int position=0;THREAD_PACKAGE pack(size);pack.send_tid=rank;pack.recv_tid=all_neighbor_ranks(n);
        for(int axis=1;axis<=T_GRID::dimension;axis++) for(FACE_ITERATOR iterator(local_grid,send_regions_n(axis),axis);iterator.Valid();iterator.Next()) 
            data.data(axis).Pack(pack.buffer,position,iterator.Face_Index());
        pthread_mutex_lock(lock);
        buffers.Append(pack);
        pthread_mutex_unlock(lock);}
    // receive
    ARRAY<ARRAY<RANGE<TV_INT> > > recv_regions(T_GRID::dimension);
    for(int axis=1;axis<=T_GRID::dimension;axis++)Find_Boundary_Regions(recv_regions(axis),Face_Sentinels(axis),true,ghost_band,true);
    pthread_barrier_wait(barr);
    for(int n=1;n<=recv_regions(1).m;n++)if(all_neighbor_ranks(n)!=-1){int index=0;
        ARRAY<RANGE<TV_INT> > recv_regions_n(T_GRID::dimension);for(int axis=1;axis<=T_GRID::dimension;axis++)recv_regions_n(axis)=recv_regions(axis)(n);
        for(int i=1;i<=buffers.m;i++) if(buffers(i).send_tid==all_neighbor_ranks(n) && buffers(i).recv_tid==rank) index=i;
        assert(index);int position=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++) for(FACE_ITERATOR iterator(local_grid,recv_regions_n(axis),axis);iterator.Valid();iterator.Next())
            data.data(axis).Unpack(buffers(index).buffer,position,iterator.Face_Index());}
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Average_Common_Face_Data
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Average_Common_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data) const
{
#ifdef USE_PTHREADS
    ARRAY<ARRAY<RANGE<TV_INT> > > regions(T_GRID::dimension);
    for(int axis=1;axis<=T_GRID::dimension;axis++)Find_Boundary_Regions(regions(axis),Face_Sentinels(axis),false,RANGE<VECTOR<int,1> >(0,0),false);
    // send and receive into temporary buffers
    for(int n=1;n<=regions(1).m;n++)if(side_neighbor_ranks(n)!=-1){int axis=(n-1)/2+1;
        int size=0;for(FACE_ITERATOR iterator(local_grid,regions(axis)(n),axis);iterator.Valid();iterator.Next()) size+=data.data(axis).Pack_Size();
        int position=0;THREAD_PACKAGE pack(size);pack.send_tid=rank;pack.recv_tid=side_neighbor_ranks(n);
        for(FACE_ITERATOR iterator(local_grid,regions(axis)(n),axis);iterator.Valid();iterator.Next()) data.data(axis).Pack(pack.buffer,position,iterator.Face_Index());
        pthread_mutex_lock(lock);
        buffers.Append(pack);
        pthread_mutex_unlock(lock);}
    // wait
    pthread_barrier_wait(barr);
    // average received data with local data (TODO: find a cleaner general way to do this)
    for(int n=1;n<=regions(1).m;n++)if(side_neighbor_ranks(n)!=-1){int axis=(n-1)/2+1;
        ARRAY<T2,FACE_INDEX<TV::dimension> > local_data;local_data.Resize(data.Domain_Indices(),false,false);
        int index=0;for(int i=1;i<=buffers.m;i++) if(buffers(i).send_tid==side_neighbor_ranks(n) && buffers(i).recv_tid==rank) index=i;
        assert(index);int position=0;
        for(FACE_ITERATOR iterator(local_grid,regions(axis)(n),axis);iterator.Valid();iterator.Next()){
            local_data.data(axis).Unpack(buffers(index).buffer,position,iterator.Face_Index());
            data(iterator.Full_Index())=data(iterator.Full_Index())*0.5+local_data(iterator.Full_Index())*0.5;}}
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Assert_Common_Face_Data
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Assert_Common_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data,const T tolerance) const
{
#ifdef USE_PTHREADS
    ARRAY<ARRAY<RANGE<TV_INT> > > regions(T_GRID::dimension);
    for(int axis=1;axis<=T_GRID::dimension;axis++)Find_Boundary_Regions(regions(axis),Face_Sentinels(axis),false,RANGE<VECTOR<int,1> >(0,0),false);
    // send and receive into temporary buffers
    for(int n=1;n<=regions(1).m;n++)if(side_neighbor_ranks(n)!=-1){int axis=(n-1)/2+1;
        int size=0;for(FACE_ITERATOR iterator(local_grid,regions(axis)(n),axis);iterator.Valid();iterator.Next()) size+=data.data(axis).Pack_Size();
        int position=0;THREAD_PACKAGE pack(size);pack.send_tid=rank;pack.recv_tid=side_neighbor_ranks(n);
        for(FACE_ITERATOR iterator(local_grid,regions(axis)(n),axis);iterator.Valid();iterator.Next()) data.data(axis).Pack(pack.buffer,position,iterator.Face_Index());
        pthread_mutex_lock(lock);
        buffers.Append(pack);
        pthread_mutex_unlock(lock);}
    // wait
    pthread_barrier_wait(barr);
    // average received data with local data (TODO: find a cleaner general way to do this)
    for(int n=1;n<=regions(1).m;n++)if(side_neighbor_ranks(n)!=-1){int axis=(n-1)/2+1;
        ARRAY<T2,FACE_INDEX<TV::dimension> > local_data;local_data.Resize(data.Domain_Indices(),false,false);
        int index=0;for(int i=1;i<=buffers.m;i++) if(buffers(i).send_tid==side_neighbor_ranks(n) && buffers(i).recv_tid==rank) index=i;
        assert(index);int position=0;
        for(FACE_ITERATOR iterator(local_grid,regions(axis)(n),axis);iterator.Valid();iterator.Next()){
            local_data.data(axis).Unpack(buffers(index).buffer,position,iterator.Face_Index());
            assert(abs(data(iterator.Full_Index())-local_data(iterator.Full_Index()))<=tolerance);}}
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
//#####################################################################
// Function Sync_Scalar
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Sync_Scalar(const ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& local_data,ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& global_data) const
{
    for(CELL_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();global_data(cell+local_to_global_offset)=local_data(cell);}
}
//#####################################################################
// Function Distribute_Scalar
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Distribute_Scalar(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& local_data,const ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& global_data) const
{
    for(CELL_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();local_data(cell)=global_data(cell+local_to_global_offset);}
}
//#####################################################################
// Function Sync_Face_Scalar
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Sync_Face_Scalar(const ARRAY<T2,FACE_INDEX<TV::dimension> >& local_data,ARRAY<T2,FACE_INDEX<TV::dimension> >& global_data) const
{
    for(int axis=1;axis<=TV::dimension;axis++){
        RANGE<TV_INT> domain=local_grid.Domain_Indices();if(domain.max_corner(axis)+local_to_global_offset(axis)==global_grid.Domain_Indices().max_corner(axis)) domain.max_corner(axis)++;
        for(FACE_ITERATOR iterator(local_grid,domain,axis);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();
            global_data.Component(axis)(face+local_to_global_offset)=local_data.Component(axis)(face);}}
}
//#####################################################################
// Function Distribute_Face_Scalar
//#####################################################################
template<class T_GRID> template<class T2> void THREADED_UNIFORM_GRID<T_GRID>::
Distribute_Face_Scalar(ARRAY<T2,FACE_INDEX<TV::dimension> >& local_data,const ARRAY<T2,FACE_INDEX<TV::dimension> >& global_data) const
{
    for(FACE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
        local_data.Component(axis)(face)=global_data.Component(axis)(face+local_to_global_offset);}
}
//#####################################################################
// Function Allgather
//#####################################################################
template<class T_GRID> void THREADED_UNIFORM_GRID<T_GRID>::
Allgather(ARRAY<int>& data) const
{
#ifdef USE_PTHREADS
    THREAD_PACKAGE pack(sizeof(int));*(int*)(&pack.buffer(1))=data(tid);pack.send_tid=tid;
    pthread_mutex_lock(lock);
    buffers.Append(pack);
    pthread_mutex_unlock(lock);
    pthread_barrier_wait(barr);
    for(int i=1;i<=buffers.m;i++) data(buffers(i).send_tid)=*(int*)(&buffers(i).buffer(1));
    pthread_barrier_wait(barr);
    if(tid==1) buffers.m=0;
    pthread_barrier_wait(barr);
#endif
}
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER(T,T_GRID,d)      \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<int>::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<bool>::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,1> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,2> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,3> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,4> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,5> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<MATRIX<T,1> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<SYMMETRIC_MATRIX<T,2> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Cell_Data(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<SYMMETRIC_MATRIX<T,3> >::TYPE&,const int,const bool) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<bool>::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,1> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,2> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,3> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,4> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,5> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<MATRIX<T,1> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,2> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Exchange_Boundary_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,3> >::TYPE&,const int) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<bool>::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,1> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,2> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,3> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,4> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<VECTOR<T,5> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<MATRIX<T,1> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,2> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Average_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,3> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Assert_Common_Face_Data(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&,const T) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Sync_Scalar(const GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE&,GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Sync_Scalar(const GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<bool>::TYPE&,GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<bool>::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Sync_Scalar(const GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,d> >::TYPE&,GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,d> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Distribute_Scalar(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE&,const GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Distribute_Scalar(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<bool>::TYPE&,const GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<bool>::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Distribute_Scalar(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,d> >::TYPE&,const GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::REBIND<VECTOR<T,d> >::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Sync_Face_Scalar(const GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&,GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Sync_Face_Scalar(const GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<bool>::TYPE&,GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<bool>::TYPE&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Distribute_Face_Scalar(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&,const GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&) const; \
    template void THREADED_UNIFORM_GRID<T_GRID>::Distribute_Face_Scalar(GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<bool>::TYPE&,const GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::REBIND<bool>::TYPE&) const; \
    template class THREADED_UNIFORM_GRID<T_GRID>;
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,1> >),1);
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,2> >),2);
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,3> >),3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,1> >),1);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,2> >),2);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,3> >),3);
#endif
