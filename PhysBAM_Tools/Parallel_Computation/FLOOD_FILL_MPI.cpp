//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/FLOOD_FILL_MPI.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FLOOD_FILL_MPI<T_GRID>::
FLOOD_FILL_MPI(const T_MPI_GRID& mpi_grid_input,const T_GRID& local_grid_input,const T_FACE_ARRAYS_BOOL& psi_N_input,int& number_of_regions_input,T_ARRAYS_INT& colors_input,
    ARRAY<ARRAY<int> >& color_ranks_input,ARRAY<bool>* color_touches_uncolorable_input)
    :mpi_grid(mpi_grid_input),local_grid(local_grid_input),psi_N(psi_N_input),number_of_regions(number_of_regions_input),colors(colors_input),color_ranks(color_ranks_input),
    color_touches_uncolorable(color_touches_uncolorable_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FLOOD_FILL_MPI<T_GRID>::
~FLOOD_FILL_MPI()
{}
//#####################################################################

#ifdef USE_MPI

//#####################################################################
// Function Union_Find_Merge_Op
//#####################################################################
static MPI::Comm* union_find_merge_op_comm; // TODO: This is needed because user-defined operations don't get to know their communicator.  A better hack would be nice.
static void Union_Find_Merge_Op(const void* in,void* inout,int len,const MPI::Datatype& datatype)
{
    // unpack union finds
    UNION_FIND<> union_find_in,union_find_inout;
    int m=datatype.Get_size();assert(len==1);
    {int position=0;MPI_UTILITIES::Unpack(union_find_in,ARRAY_VIEW<const char>(m,static_cast<const char*>(in)),position,*union_find_merge_op_comm);}
    {int position=0;MPI_UTILITIES::Unpack(union_find_inout,ARRAY_VIEW<const char>(m,static_cast<char*>(inout)),position,*union_find_merge_op_comm);}
    // merge
    union_find_inout.Merge(union_find_in);
    // pack output 
    {int position=0;MPI_UTILITIES::Pack(union_find_inout,ARRAY_VIEW<char>(m,static_cast<char*>(inout)),position,*union_find_merge_op_comm);}
}
//#####################################################################
// Function Synchronize_Colors
//#####################################################################
template<class T_ARRAYS,class T_GRID,class T_BOX> static inline void Resize_Helper(T_ARRAYS& array,const T_GRID& grid,const T_BOX& box)
{
    array.Resize(box); // for uniform grids we can allocate exactly the requested region
}
template<class T_GRID,class T_BOX> static inline void Resize_Helper(ARRAY<int>& array,const T_GRID& grid,const T_BOX& box)
{
    array.Resize(grid.Cell_Indices()); // for rle grids we have to allocate space for the entire grid due to the unstructuredness
}
template<class T_GRID> int FLOOD_FILL_MPI<T_GRID>::
Synchronize_Colors()
{
    if(mpi_grid.threaded_grid) return Synchronize_Colors_Threaded();
    ARRAY<RANGE<typename T_PARALLEL_GRID::VECTOR_INT> > boundary_regions;
    mpi_grid.Find_Boundary_Regions(boundary_regions,RANGE<typename T_PARALLEL_GRID::VECTOR_INT>::Zero_Box(),false,RANGE<VECTOR<int,1> >(-1,0),false,true,local_grid);
    // figure out which colors are global
    int global_color_count=0;
    ARRAY<int,VECTOR<int,1> > color_map(-1,number_of_regions);color_map(-1)=-1;color_map(0)=0;
    {ARRAY<bool,VECTOR<int,1> > color_is_global(-1,number_of_regions);
    Find_Global_Colors(color_is_global,RANGE<typename T_PARALLEL_GRID::VECTOR_INT>::Centered_Box());
    for(int color=1;color<=number_of_regions;color++)if(color_is_global(color)) color_map(color)=++global_color_count;}

    // send numbers of global colors to everyone
    ARRAY<int> global_color_counts(mpi_grid.number_of_processes);
    mpi_grid.comm->Allgather(&global_color_count,1,MPI_UTILITIES::Datatype<int>(),&global_color_counts(1),1,MPI_UTILITIES::Datatype<int>());
    int total_global_colors=ARRAYS_COMPUTATIONS::Sum(global_color_counts);
    int global_color_offset=ARRAYS_COMPUTATIONS::Sum(global_color_counts.Prefix(mpi_grid.rank));
    LOG::cout<<"initial colors: "<<number_of_regions<<" total, "<<global_color_count<<" out of "<<total_global_colors<<" global"<<std::endl;
    if(!total_global_colors){color_ranks.Clean_Memory();return 0;}

    ARRAY<MPI_PACKAGE> packages;
    ARRAY<T_ARRAYS_INT> colors_copy(boundary_regions.m);
    // send left (front) colors
    ARRAY<MPI::Request> send_requests;
    for(int side=1;side<=T_PARALLEL_GRID::number_of_faces_per_cell;side+=2)if(mpi_grid.side_neighbor_ranks(side)!=MPI::PROC_NULL){
        Resize_Helper(colors_copy(side),local_grid,boundary_regions(side));
        Translate_Local_Colors_To_Global_Colors(color_map,colors_copy(side),boundary_regions(side),global_color_offset);
        MPI_PACKAGE package=mpi_grid.Package_Cell_Data(colors_copy(side),boundary_regions(side));
        packages.Append(package);
        send_requests.Append(package.Isend(*mpi_grid.comm,mpi_grid.side_neighbor_ranks(side),mpi_grid.Get_Send_Tag(mpi_grid.side_neighbor_directions(side))));}
    // receive right (back) colors and initialize union find
    UNION_FIND<> union_find(total_global_colors);
    {ARRAY<MPI::Request> recv_requests;
    for(int side=2;side<=T_PARALLEL_GRID::number_of_faces_per_cell;side+=2)if(mpi_grid.side_neighbor_ranks(side)!=MPI::PROC_NULL){
        Resize_Helper(colors_copy(side),local_grid,boundary_regions(side));
        MPI_PACKAGE package=mpi_grid.Package_Cell_Data(colors_copy(side),boundary_regions(side));
        packages.Append(package);
        recv_requests.Append(package.Irecv(*mpi_grid.comm,mpi_grid.side_neighbor_ranks(side),mpi_grid.Get_Recv_Tag(mpi_grid.side_neighbor_directions(side))));}
    MPI::Status status;
    while(MPI_UTILITIES::Wait_Any(recv_requests,status)){
        int side;for(side=2;side<=T_PARALLEL_GRID::number_of_faces_per_cell;side+=2)if(mpi_grid.Get_Recv_Tag(mpi_grid.side_neighbor_directions(side))==status.Get_tag()) break;
        Find_Color_Matches(color_map,union_find,colors_copy(side),boundary_regions(side),global_color_offset);}}

    // synchronize union find
    UNION_FIND<> final_union_find;
    {ARRAY<char> union_find_buffer(MPI_UTILITIES::Pack_Size(union_find,*mpi_grid.comm)+1);
    {int position=0;MPI_UTILITIES::Pack(union_find,union_find_buffer,position,*mpi_grid.comm);}
    MPI::Datatype union_find_type=MPI::PACKED.Create_contiguous(union_find_buffer.m);union_find_type.Commit();
    MPI::Op union_find_merge_op;union_find_merge_op.Init(Union_Find_Merge_Op,true);
    ARRAY<char> final_union_find_buffer(union_find_buffer.m);
    union_find_merge_op_comm=mpi_grid.comm;
    mpi_grid.comm->Allreduce(union_find_buffer.Get_Array_Pointer(),final_union_find_buffer.Get_Array_Pointer(),1,union_find_type,union_find_merge_op);
    {int position=0;MPI_UTILITIES::Unpack(final_union_find,final_union_find_buffer,position,*mpi_grid.comm);}
    union_find_type.Free();union_find_merge_op.Free();}

    // fix color map for global colors
    number_of_regions=0;
    ARRAY<int> global_to_final_color_map(total_global_colors);
    for(int i=1;i<=total_global_colors;i++){
        int root=final_union_find.Find(i);
        if(!global_to_final_color_map(root)) global_to_final_color_map(root)=++number_of_regions;
        global_to_final_color_map(i)=global_to_final_color_map(root);}
    for(int i=1;i<=color_map.domain.max_corner.x;i++)if(color_map(i)>0) color_map(i)=global_to_final_color_map(color_map(i)+global_color_offset);

    // find list of processes corresponding to each color
    int end=0;
    color_ranks.Clean_Memory();
    color_ranks.Resize(number_of_regions);
    for(int r=0;r<mpi_grid.number_of_processes;r++){
        int start=end+1;end+=global_color_counts(r+1);
        for(int i=start;i<=end;i++)color_ranks(global_to_final_color_map(i)).Append_Unique(r);}
    for(int color=1;color<=color_ranks.m;color++) assert(color_ranks(color).m>1 || mpi_grid.side_neighbor_ranks.Contains(mpi_grid.rank));

    // remap colors
    Remap_Colors(color_map,RANGE<typename T_PARALLEL_GRID::VECTOR_INT>::Centered_Box());

    LOG::cout<<"final colors: "<<color_ranks.m<<" global, "<<number_of_regions-color_ranks.m<<" local"<<std::endl;

    // remap color_touches_uncolorable
    if(color_touches_uncolorable){
        ARRAY<bool> new_color_touches_uncolorable(number_of_regions);
        for(int i=1;i<=color_touches_uncolorable->m;i++)if(color_map(i)>0) new_color_touches_uncolorable(color_map(i))|=(*color_touches_uncolorable)(i);
        color_touches_uncolorable->Exchange(new_color_touches_uncolorable);
        // synchronize color_touches_uncolorable, TODO: this could be merged with above communication
        ARRAY<bool> global_color_touches_uncolorable(color_ranks.m);
        ARRAY<bool>::Get(global_color_touches_uncolorable,*color_touches_uncolorable);
        mpi_grid.comm->Allreduce(&global_color_touches_uncolorable(1),&(*color_touches_uncolorable)(1),color_ranks.m,MPI_UTILITIES::Datatype<bool>(),MPI::LOR);}

    // finish
    MPI_UTILITIES::Wait_All(send_requests);
    MPI_PACKAGE::Free_All(packages);

    return color_ranks.m;
}
template<class T_GRID> int FLOOD_FILL_MPI<T_GRID>::
Synchronize_Colors_Threaded()
{
#ifdef USE_PTHREADS
    ARRAY<RANGE<typename T_PARALLEL_GRID::VECTOR_INT> > boundary_regions;
    mpi_grid.threaded_grid->Find_Boundary_Regions(boundary_regions,RANGE<typename T_PARALLEL_GRID::VECTOR_INT>::Zero_Box(),false,RANGE<VECTOR<int,1> >(-1,0),false,true,local_grid);
    // figure out which colors are global
    int global_color_count=0;
    ARRAY<int,VECTOR<int,1> > color_map(-1,number_of_regions);color_map(-1)=-1;color_map(0)=0;
    {ARRAY<bool,VECTOR<int,1> > color_is_global(-1,number_of_regions);
    Find_Global_Colors(color_is_global,RANGE<typename T_PARALLEL_GRID::VECTOR_INT>::Centered_Box());
    for(int color=1;color<=number_of_regions;color++)if(color_is_global(color)) color_map(color)=++global_color_count;}

    // send numbers of global colors to everyone
    ARRAY<int> global_color_counts(mpi_grid.threaded_grid->number_of_processes);
    mpi_grid.threaded_grid->Allgather(global_color_counts);
    int total_global_colors=ARRAYS_COMPUTATIONS::Sum(global_color_counts);
    int global_color_offset=ARRAYS_COMPUTATIONS::Sum(global_color_counts.Prefix(mpi_grid.threaded_grid->rank));
    LOG::cout<<"initial colors: "<<number_of_regions<<" total, "<<global_color_count<<" out of "<<total_global_colors<<" global"<<std::endl;
    if(!total_global_colors){color_ranks.Clean_Memory();return 0;}

    ARRAY<T_ARRAYS_INT> colors_copy(boundary_regions.m);
    // send left (front) colors
    for(int side=1;side<=T_PARALLEL_GRID::number_of_faces_per_cell;side+=2)if(mpi_grid.threaded_grid->side_neighbor_ranks(side)!=-1){
        Resize_Helper(colors_copy(side),local_grid,boundary_regions(side));
        Translate_Local_Colors_To_Global_Colors(color_map,colors_copy(side),boundary_regions(side),global_color_offset);
        THREAD_PACKAGE pack=mpi_grid.threaded_grid->Package_Cell_Data(colors_copy(side),boundary_regions(side));pack.recv_tid=mpi_grid.threaded_grid->side_neighbor_ranks(side);
        pthread_mutex_lock(mpi_grid.threaded_grid->lock);
        mpi_grid.threaded_grid->buffers.Append(pack);
        pthread_mutex_unlock(mpi_grid.threaded_grid->lock);}
    // receive right (back) colors and initialize union find
    UNION_FIND<> union_find(total_global_colors);
    pthread_barrier_wait(mpi_grid.threaded_grid->barr);
    {for(int side=2;side<=T_PARALLEL_GRID::number_of_faces_per_cell;side+=2)if(mpi_grid.threaded_grid->side_neighbor_ranks(side)!=-1){
        Resize_Helper(colors_copy(side),local_grid,boundary_regions(side));
        int index=0;for(int i=1;i<=mpi_grid.threaded_grid->buffers.m;i++) if(mpi_grid.threaded_grid->buffers(i).send_tid==mpi_grid.threaded_grid->side_neighbor_ranks(side) && mpi_grid.threaded_grid->buffers(i).recv_tid==mpi_grid.threaded_grid->rank) index=i;
        PHYSBAM_ASSERT(index);int position=0;
        for(CELL_ITERATOR iterator(local_grid,boundary_regions(side));iterator.Valid();iterator.Next()) colors_copy(side).Unpack(mpi_grid.threaded_grid->buffers(index).buffer,position,iterator.Cell_Index());
        Find_Color_Matches(color_map,union_find,colors_copy(side),boundary_regions(side),global_color_offset);}}
    pthread_barrier_wait(mpi_grid.threaded_grid->barr);
    if(mpi_grid.threaded_grid->tid==1) mpi_grid.threaded_grid->buffers.m=0;
    pthread_barrier_wait(mpi_grid.threaded_grid->barr);

    // synchronize union find
    UNION_FIND<> final_union_find;
    {THREAD_PACKAGE pack(sizeof(int)*(2+union_find.parents.m)+union_find.ranks.m);
    {int position=0;*(int*)(&pack.buffer(position+1))=union_find.parents.m;position+=sizeof(int);*(int*)(&pack.buffer(position+1))=union_find.ranks.m;position+=sizeof(int);
    for(int i=1;i<=union_find.parents.m;i++) union_find.parents.Pack(pack.buffer,position,i);for(int i=1;i<=union_find.ranks.m;i++) union_find.ranks.Pack(pack.buffer,position,i);
    pthread_mutex_lock(mpi_grid.threaded_grid->lock);
    mpi_grid.threaded_grid->buffers.Append(pack);
    pthread_mutex_unlock(mpi_grid.threaded_grid->lock);}
    pthread_barrier_wait(mpi_grid.threaded_grid->barr);
    for(int buf=1;buf<=mpi_grid.threaded_grid->buffers.m;buf++){int position=0;
        union_find.parents.Resize(*(int*)(&pack.buffer(position+1)));position+=sizeof(int);union_find.ranks.Resize(*(int*)(&pack.buffer(position+1)));position+=sizeof(int);
        for(int i=1;i<=union_find.parents.m;i++) union_find.parents.Unpack(mpi_grid.threaded_grid->buffers(buf).buffer,position,i);for(int i=1;i<=union_find.ranks.m;i++) union_find.ranks.Unpack(mpi_grid.threaded_grid->buffers(buf).buffer,position,i);
        final_union_find.Merge(union_find);}
    pthread_barrier_wait(mpi_grid.threaded_grid->barr);
    if(mpi_grid.threaded_grid->tid==1) mpi_grid.threaded_grid->buffers.m=0;
    pthread_barrier_wait(mpi_grid.threaded_grid->barr);}

    // fix color map for global colors
    number_of_regions=0;
    ARRAY<int> global_to_final_color_map(total_global_colors);
    for(int i=1;i<=total_global_colors;i++){
        int root=final_union_find.Find(i);
        if(!global_to_final_color_map(root)) global_to_final_color_map(root)=++number_of_regions;
        global_to_final_color_map(i)=global_to_final_color_map(root);}
    for(int i=1;i<=color_map.domain.max_corner.x;i++)if(color_map(i)>0) color_map(i)=global_to_final_color_map(color_map(i)+global_color_offset);

    // find list of processes corresponding to each color
    int end=0;
    color_ranks.Clean_Memory();
    color_ranks.Resize(number_of_regions);
    for(int r=0;r<mpi_grid.number_of_processes;r++){
        int start=end+1;end+=global_color_counts(r+1);
        for(int i=start;i<=end;i++)color_ranks(global_to_final_color_map(i)).Append_Unique(r);}
    for(int color=1;color<=color_ranks.m;color++) assert(color_ranks(color).m>1 || mpi_grid.side_neighbor_ranks.Contains(mpi_grid.rank));

    // remap colors
    Remap_Colors(color_map,RANGE<typename T_PARALLEL_GRID::VECTOR_INT>::Centered_Box());

    LOG::cout<<"final colors: "<<color_ranks.m<<" global, "<<number_of_regions-color_ranks.m<<" local"<<std::endl;

    // remap color_touches_uncolorable
    if(color_touches_uncolorable){
        ARRAY<bool> new_color_touches_uncolorable(number_of_regions);
        for(int i=1;i<=color_touches_uncolorable->m;i++)if(color_map(i)>0) new_color_touches_uncolorable(color_map(i))|=(*color_touches_uncolorable)(i);
        color_touches_uncolorable->Exchange(new_color_touches_uncolorable);
        // synchronize color_touches_uncolorable, TODO: this could be merged with above communication
        ARRAY<bool> global_color_touches_uncolorable(color_ranks.m);
        ARRAY<bool>::Get(global_color_touches_uncolorable,*color_touches_uncolorable);
        THREAD_PACKAGE pack(sizeof(bool)*global_color_touches_uncolorable.m);
        {int position=0;for(int i=1;i<=global_color_touches_uncolorable.m;i++) global_color_touches_uncolorable.Pack(pack.buffer,position,i);
        pthread_mutex_lock(mpi_grid.threaded_grid->lock);
        mpi_grid.threaded_grid->buffers.Append(pack);
        pthread_mutex_unlock(mpi_grid.threaded_grid->lock);}
        pthread_barrier_wait(mpi_grid.threaded_grid->barr);
        for(int buf=1;buf<=mpi_grid.threaded_grid->buffers.m;buf++){int position=0;
            for(int i=1;i<=global_color_touches_uncolorable.m;i++) global_color_touches_uncolorable.Unpack(mpi_grid.threaded_grid->buffers(buf).buffer,position,i);
            for(int i=1;i<=global_color_touches_uncolorable.m;i++) (*color_touches_uncolorable)(i)|=global_color_touches_uncolorable(i);}
        pthread_barrier_wait(mpi_grid.threaded_grid->barr);
        if(mpi_grid.threaded_grid->tid==1) mpi_grid.threaded_grid->buffers.m=0;
        pthread_barrier_wait(mpi_grid.threaded_grid->barr);}

    return color_ranks.m;
#endif
}

//#####################################################################
// Function Find_Global_Colors
//#####################################################################
template<class T_GRID> void FLOOD_FILL_MPI<T_GRID>::
Find_Global_Colors(ARRAY<bool,VECTOR<int,1> >& color_is_global,const RANGE<TV_INT>&) const
{
    int proc_null=mpi_grid.threaded_grid?-1:MPI::PROC_NULL;
    const ARRAY<int>& side_neighbor_ranks=mpi_grid.threaded_grid?mpi_grid.threaded_grid->side_neighbor_ranks:mpi_grid.side_neighbor_ranks;
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=0;axis_side<=1;axis_side++){
        int side=2*axis-1+axis_side;
        if(side_neighbor_ranks(side)!=proc_null){
            for(typename T_GRID::FACE_ITERATOR iterator(local_grid,0,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next())if(!psi_N.Component(axis)(iterator.Face_Index())){
                color_is_global(colors(iterator.First_Cell_Index()))=true;
                color_is_global(colors(iterator.Second_Cell_Index()))=true;}}}
}
//#####################################################################
// Function Find_Global_Colors
//#####################################################################
template<class T_GRID> template<class T_BOX_HORIZONTAL_INT> void FLOOD_FILL_MPI<T_GRID>::
Find_Global_Colors(ARRAY<bool,VECTOR<int,1> >& color_is_global,const T_BOX_HORIZONTAL_INT&) const
{
    T_GRID::template Horizontal_Face_Loop<Find_Global_Colors_Helper>(mpi_grid,mpi_grid.local_grid,psi_N,colors,color_is_global);
}
//#####################################################################
// Function Find_Global_Colors_Helper
//#####################################################################
template<class T_GRID> template<class T_FACE> void FLOOD_FILL_MPI<T_GRID>::
Find_Global_Colors_Helper::Apply(const T_MPI_GRID& mpi_grid,const T_GRID& local_grid,const ARRAY<bool>& psi_N,const ARRAY<int>& colors,ARRAY<bool,VECTOR<int,1> >& color_is_global)
{
    int horizontal_axis=T_FACE::Horizontal_Axis();
    ARRAY<typename T_GRID::BOX_HORIZONTAL_INT> boundary_regions;
    mpi_grid.Find_Boundary_Regions(boundary_regions,CELL_ITERATOR::Sentinels(),false,RANGE<VECTOR<int,1> >(0,0),false,true,local_grid);
    for(int axis_side=0;axis_side<=1;axis_side++){
        int side=2*horizontal_axis-1+axis_side;
        if(mpi_grid.side_neighbor_ranks(side)!=MPI::PROC_NULL){
            typename T_GRID::BOX_HORIZONTAL_INT face_region=boundary_regions(side);if(axis_side) face_region+=T_GRID::VECTOR_HORIZONTAL_INT::Axis_Vector(horizontal_axis);
            for(T_FACE face(local_grid,face_region);face;face++)if(!psi_N(face.Face())){int c1=face.cell1.Cell(),c2=face.cell2.Cell();
                for(int i=0;i<=face.cell1.Long();i++)color_is_global(colors(c1+i))=true;
                for(int i=0;i<=face.cell2.Long();i++)color_is_global(colors(c2+i))=true;}}}
}
//#####################################################################
// Function Translate_Local_Colors_To_Global_Colors
//#####################################################################
template<class T_GRID> void FLOOD_FILL_MPI<T_GRID>::
Translate_Local_Colors_To_Global_Colors(const ARRAY<int,VECTOR<int,1> >& color_map,T_ARRAYS_INT& colors_copy,const RANGE<TV_INT>& region,const int global_color_offset) const
{
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        int new_color=color_map(colors(cell_index));colors_copy(cell_index)=new_color>0?new_color+global_color_offset:-1;}
}
//#####################################################################
// Function Translate_Local_Colors_To_Global_Colors
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T_GRID> template<class T_BOX_HORIZONTAL_INT> void FLOOD_FILL_MPI<T_GRID>::
Translate_Local_Colors_To_Global_Colors(const ARRAY<int,VECTOR<int,1> >& color_map,ARRAY<int>& colors_copy,const T_BOX_HORIZONTAL_INT& region,const int global_color_offset) const
{
    for(RLE_GRID_SIMPLE_ITERATOR<T_GRID,CELL_ITERATOR> cell(local_grid,region);cell;cell++){int c=cell.index;
        int new_color=color_map(colors(c));colors_copy(c)=new_color>0?new_color+global_color_offset:-1;}
}
#endif
//#####################################################################
// Function Find_Color_Matches
//#####################################################################
template<class T_GRID> void FLOOD_FILL_MPI<T_GRID>::
Find_Color_Matches(const ARRAY<int,VECTOR<int,1> >& color_map,UNION_FIND<>& union_find,T_ARRAYS_INT& colors_copy,const RANGE<TV_INT>& region,const int global_color_offset) const
{
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();if(colors_copy(cell_index)>0){
        int local_color=color_map(colors(cell_index));
        if(local_color>0) union_find.Union(global_color_offset+local_color,colors_copy(cell_index));}}
}
//#####################################################################
// Function Find_Color_Matches
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T_GRID> template<class T_BOX_HORIZONTAL_INT> void FLOOD_FILL_MPI<T_GRID>::
Find_Color_Matches(const ARRAY<int,VECTOR<int,1> >& color_map,UNION_FIND<>& union_find,ARRAY<int>& colors_copy,const T_BOX_HORIZONTAL_INT& region,const int global_color_offset) const
{
    for(RLE_GRID_SIMPLE_ITERATOR<T_GRID,CELL_ITERATOR> cell(local_grid,region);cell;cell++){int color=colors_copy(cell.index);if(color>0){
        int local_color=color_map(colors(cell.index));
        if(local_color>0) union_find.Union(global_color_offset+local_color,color);}}
}
#endif
//#####################################################################
// Function Remap_Colors
//#####################################################################
template<class T_GRID> void FLOOD_FILL_MPI<T_GRID>::
Remap_Colors(ARRAY<int,VECTOR<int,1> >& color_map,const RANGE<TV_INT>&)
{
    for(CELL_ITERATOR iterator(local_grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        int color=colors(cell_index);assert(color); // colors should either be -1 or positive
        if(color_map(color)==0) color_map(color)=++number_of_regions;
        colors(cell_index)=color_map(color);}
}
//#####################################################################
// Function Remap_Colors
//#####################################################################
template<class T_GRID> template<class T_BOX_HORIZONTAL_INT> void FLOOD_FILL_MPI<T_GRID>::
Remap_Colors(ARRAY<int,VECTOR<int,1> >& color_map,const T_BOX_HORIZONTAL_INT&)
{
    for(int c=1;c<=colors.m;c++){
        int color=colors(c);assert(color); // colors should either be -1 or positive
        if(color_map(color)==0) color_map(color)=++number_of_regions;
        colors(c)=color_map(color);}
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID> int FLOOD_FILL_MPI<T_GRID>::Synchronize_Colors(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();return 0;}
//#####################################################################

#endif

//#####################################################################
#define INSTANTIATION_HELPER(T,T_GRID) \
    template FLOOD_FILL_MPI<T_GRID >::FLOOD_FILL_MPI(const T_MPI_GRID&,const T_GRID&,const T_FACE_ARRAYS_BOOL&,int&,T_ARRAYS_INT&,ARRAY<ARRAY<int> >&,ARRAY<bool>*); \
    template FLOOD_FILL_MPI<T_GRID >::~FLOOD_FILL_MPI(); \
    template int FLOOD_FILL_MPI<T_GRID >::Synchronize_Colors();
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,1> >));
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,2> >));
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,3> >));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,1> >));
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,2> >));
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,3> >));
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
INSTANTIATION_HELPER(float,RLE_GRID_2D<float>);
INSTANTIATION_HELPER(float,RLE_GRID_3D<float>);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,RLE_GRID_2D<double>);
INSTANTIATION_HELPER(double,RLE_GRID_3D<double>);
#endif
#endif
