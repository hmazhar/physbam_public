//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
using namespace PhysBAM;
//#####################################################################
// Function Serial_Solve
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Serial_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const int tag,const T tolerance)
{
    // TODO: this routine is useful only for testing purposes
    LOG::SCOPE scope("MPI SOLVE","mpi solve");
    int processors=comm.Get_size(),rank=comm.Get_rank(),master=0;assert(processors>1);
    if(rank!=master){
        // send linear system piece
        int buffer_size=MPI_UTILITIES::Pack_Size(partition,A,x,b,comm);
        ARRAY<char> buffer(buffer_size);int position=0;
        MPI_UTILITIES::Pack(partition,A,x,b,buffer,position,comm);
        comm.Send(&buffer(1),position,MPI::PACKED,master,tag);
        // receive result
        comm.Recv(&x(1),x.n,MPI_UTILITIES::Datatype<T>(),master,tag);}
    else{
        // receive linear system pieces
        ARRAY<SPARSE_MATRIX_PARTITION> partition_array(processors);partition_array(1)=partition; // TODO: very inefficient to copy everything into one array
        ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array(processors);A_array(1)=A;
        ARRAY<VECTOR_ND<T> > x_array(processors),b_array(processors);x_array(1)=x;b_array(1)=b;
        for(int p=2;p<=processors;p++){
            MPI::Status status;
            comm.Probe(MPI::ANY_SOURCE,tag,status);
            int source=status.Get_source();
            ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
            comm.Recv(&buffer(1),buffer.m,MPI::PACKED,source,tag);
            MPI_UTILITIES::Unpack(partition_array(source+1),A_array(source+1),x_array(source+1),b_array(source+1),buffer,position,comm);}
        // find total sizes and offsets
        int global_rows=0,global_entries=0;
        for(int p=1;p<=processors;p++){
            partition_array(p).Set_Interior_Offset(global_rows);
            global_rows+=partition_array(p).Interior_Rows();
            global_entries+=partition_array(p).Interior_Entries(A_array(p));}
        for(int p=1;p<=processors;p++){SPARSE_MATRIX_PARTITION& partition=partition_array(p);
            for(int s=1;s<=partition.number_of_sides;s++)
                if(partition.neighbor_ranks(s)!=MPI::PROC_NULL) partition.neighbors(s)=&partition_array(partition.neighbor_ranks(s)+1);}
        // find global offsets
        SPARSE_MATRIX_FLAT_NXN<T> global_A;
        global_A.n=global_rows;global_A.offsets.Resize(global_rows+1);
        {int current_row=1,current_index=1;global_A.offsets(1)=1;
        for(int p=1;p<=processors;p++) for(int i=partition_array(p).interior_indices.min_corner;i<=partition_array(p).interior_indices.max_corner;i++)
            global_A.offsets(++current_row)=current_index+=A_array(p).offsets(i+1)-A_array(p).offsets(i);
        assert(current_row==global_rows+1 && global_A.offsets(current_row)==global_entries+1);}
        // assemble full linear system
        global_A.A.Resize(global_entries);
        VECTOR_ND<T> global_x(global_rows),global_b(global_rows);
        for(int p=1;p<=processors;p++) for(int i=partition_array(p).interior_indices.min_corner;i<=partition_array(p).interior_indices.max_corner;i++){
            int global_i=partition_array(p).Translate_Index(i);
            global_x(global_i)=x_array(p)(i);
            global_b(global_i)=b_array(p)(i);
            for(int index=A_array(p).offsets(i);index<A_array(p).offsets(i+1);index++)
                global_A(global_i,partition_array(p).Translate_Index(A_array(p).A(index).j))=A_array(p).A(index).a;}
        // solve linear system
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(pcg.show_results) LOG::cout<<"solving "<<global_A.n<<" cells to tolerance "<<tolerance<<std::endl;
#endif
        assert(global_A.Symmetric());global_A.Positive_Diagonal_And_Nonnegative_Row_Sum((T)1e-4);
        pcg.Solve(global_A,global_x,global_b,q,s,r,k,z,tolerance);
        // take apart result
        for(int p=1;p<=processors;p++){
            SPARSE_MATRIX_PARTITION& partition=partition_array(p);
            for(int i=partition.interior_indices.min_corner;i<=partition.interior_indices.max_corner;i++) x_array(p)(i)=global_x(partition.Translate_Interior_Index(i));
            for(int region=1;region<=partition.number_of_sides;region++){
                for(int i=partition.ghost_indices(region).min_corner;i<=partition.ghost_indices(region).max_corner;i++) x_array(p)(i)=global_x(partition.Translate_Ghost_Index(i,region));}}
        // send result pieces
        ARRAY<MPI::Request> requests;
        for(int p=2;p<=processors;p++)requests.Append(comm.Isend(&x_array(p)(1),x_array(p).n,MPI_UTILITIES::Datatype<T>(),p-1,tag));
        MPI_UTILITIES::Wait_All(requests);
        x=x_array(1);}
}
//#####################################################################
// Function Parallel_Solve
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance,const bool recompute_preconditioner)
{
    if(thread_grid){Parallel_Solve(A,x,b,thread_grid->global_column_index_boundaries,tolerance,recompute_preconditioner);return;}
    Initialize_Datatypes();
    int local_n=A.n,interior_n=partition.interior_indices.Size()+1;
    int global_n=Global_Sum(interior_n);
    T global_tolerance=Global_Max(tolerance);
    int desired_iterations=global_n;if(pcg.enforce_compatibility) desired_iterations--;if(pcg.maximum_iterations) desired_iterations=min(desired_iterations,pcg.maximum_iterations);

    VECTOR_ND<T> temp(local_n,false),p(local_n,false),z_interior(interior_n,false);

    // build interior views of x,b,p,z,temp
    VECTOR_ND<T> x_interior,b_interior,p_interior,temp_interior;
    x_interior.Set_Subvector_View(x,partition.interior_indices);
    b_interior.Set_Subvector_View(b,partition.interior_indices);
    p_interior.Set_Subvector_View(p,partition.interior_indices);
    temp_interior.Set_Subvector_View(temp,partition.interior_indices);

    // adjust x for the null space
    if(pcg.enforce_compatibility && pcg.remove_null_space_solution_component) x_interior-=(T)(Global_Sum(x_interior.Sum_Double_Precision())/global_n);

    // find initial residual, r=b-Ax - reusing b for the residual
    Fill_Ghost_Cells(x);
    A.Times(x,temp);b_interior-=temp_interior;
    if(pcg.enforce_compatibility) b_interior-=(T)(Global_Sum(b_interior.Sum_Double_Precision())/global_n);
    if(Global_Max(b_interior.Max_Abs())<=global_tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(pcg.show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
#endif
        return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    if(pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)){
        delete A.C;A.C=A.Create_Submatrix(partition.interior_indices);
        A.C->In_Place_Incomplete_Cholesky_Factorization(pcg.modified_incomplete_cholesky,pcg.modified_incomplete_cholesky_coefficient,
            pcg.preconditioner_zero_tolerance,pcg.preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        if(pcg.incomplete_cholesky){
            // solve Mz=r
            A.C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
            A.C->Solve_Backward_Substitution(temp_interior,z_interior,false,true);} // diagonal is inverted to save on divides
        else z_interior=b_interior; // set z=r when there is no preconditioner

        // for Neumann boundary conditions only, make sure z sums to zero
        if(pcg.enforce_compatibility) z_interior-=(T)(Global_Sum(z_interior.Sum_Double_Precision())/global_n);

        // update search direction
        rho_old=rho;rho=Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior,b_interior));
        T beta=0;if(iteration==1) p_interior=z_interior;else{beta=(T)(rho/rho_old);for(int i=1;i<=interior_n;i++) p_interior(i)=z_interior(i)+beta*p_interior(i);} // when iteration=1, beta=0

        // update solution and residual
        Fill_Ghost_Cells(p);
        A.Times(p,temp);
        T alpha=(T)(rho/Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior)));
        for(int i=1;i<=interior_n;i++){x_interior(i)+=alpha*p_interior(i);b_interior(i)-=alpha*temp_interior(i);}

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(pcg.enforce_compatibility) b_interior-=(T)(Global_Sum(b_interior.Sum_Double_Precision())/global_n);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        T residual=Global_Max(b_interior.Max_Abs());

        // check for convergence
        if(pcg.show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(pcg.show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(pcg.show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
#endif
    }

    Fill_Ghost_Cells(x);
}
//#####################################################################
// Function Parallel_Solve
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x_local,VECTOR_ND<T>& b_local,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries,
    const T tolerance,const bool recompute_preconditioner)
{
    // TODO templatize this with the one above
    // Set up arrays to exchange with all other fluid procs
    Find_Ghost_Regions(A,proc_column_index_boundaries);
    int global_n=A.n;
    T global_tolerance=Global_Max(tolerance);
    int desired_iterations=global_n;if(pcg.enforce_compatibility) desired_iterations--;if(pcg.maximum_iterations) desired_iterations=min(desired_iterations,pcg.maximum_iterations);
    int my_rank=comm.Get_rank();

    VECTOR_ND<T> temp(global_n,false),z(global_n,false),p_interior,x_interior;
    VECTOR_ND<T> x_global(global_n),p_global(global_n,false);
    for(int i=1;i<=global_n;i++) x_global(i+proc_column_index_boundaries(my_rank+1).x-1)=x_local(i);
    p_interior.Set_Subvector_View(p_global,INTERVAL<int>(proc_column_index_boundaries(my_rank+1).x,proc_column_index_boundaries(my_rank+1).y));
    x_interior.Set_Subvector_View(x_global,INTERVAL<int>(proc_column_index_boundaries(my_rank+1).x,proc_column_index_boundaries(my_rank+1).y));

    // adjust x for the null space
    if(pcg.enforce_compatibility && pcg.remove_null_space_solution_component) x_global-=(T)(Global_Sum(x_interior.Sum_Double_Precision())/global_n);

    // find initial residual, r=b-Ax - reusing b for the residual
    Fill_Ghost_Cells_Far(x_global);
    A.Times(x_global,temp);b_local-=temp;
    if(pcg.enforce_compatibility) b_local-=(T)(Global_Sum(b_local.Sum_Double_Precision())/global_n);
    if(Global_Max(b_local.Max_Abs())<=global_tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(pcg.show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
#endif
        return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    if(pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)){
        if(A.C)delete A.C;A.C=A.Create_Submatrix(INTERVAL<int>(proc_column_index_boundaries(my_rank+1).x,proc_column_index_boundaries(my_rank+1).y));
        A.C->In_Place_Incomplete_Cholesky_Factorization(pcg.modified_incomplete_cholesky,pcg.modified_incomplete_cholesky_coefficient,
            pcg.preconditioner_zero_tolerance,pcg.preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        if(pcg.incomplete_cholesky){
            // solve Mz=r
            A.C->Solve_Forward_Substitution(b_local,temp,true); // diagonal should be treated as the identity
            A.C->Solve_Backward_Substitution(temp,z,false,true);} // diagonal is inverted to save on divides
        else z=b_local; // set z=r when there is no preconditioner

        // for Neumann boundary conditions only, make sure z sums to zero
        if(pcg.enforce_compatibility) z-=(T)(Global_Sum(z.Sum_Double_Precision())/global_n);

        // update search direction
        rho_old=rho;rho=Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(z,b_local));
        T beta=0;if(iteration==1) p_interior=z;else{beta=(T)(rho/rho_old);for(int i=1;i<=global_n;i++) p_interior(i)=z(i)+beta*p_interior(i);} // when iteration=1, beta=0

        // update solution and residual
        Fill_Ghost_Cells_Far(p_global);
        A.Times(p_global,temp);
        T alpha=(T)(rho/Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp)));
        for(int i=1;i<=global_n;i++){x_interior(i)+=alpha*p_interior(i);b_local(i)-=alpha*temp(i);}

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(pcg.enforce_compatibility) b_local-=(T)(Global_Sum(b_local.Sum_Double_Precision())/global_n);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        T residual=Global_Max(b_local.Max_Abs());

        // check for convergence
        if(pcg.show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(pcg.show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(pcg.show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
#endif
    }

    // Copy stuff back into x_local
    for(int i=1;i<=global_n;i++) x_local(i)=x_global(i+proc_column_index_boundaries(my_rank+1).x-1);
}
//#####################################################################
// Function Initialize_Datatypes
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Initialize_Datatypes()
{
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);
    boundary_datatypes.Resize(partition.number_of_sides);ghost_datatypes.Resize(partition.number_of_sides);
    for(int s=1;s<=partition.number_of_sides;s++) if(partition.neighbor_ranks(s)!=MPI::PROC_NULL){
        if(partition.boundary_indices(s).m){
            const ARRAY<int>& displacements=partition.boundary_indices(s);
            ARRAY<int> block_lengths(displacements.m,false);ARRAYS_COMPUTATIONS::Fill(block_lengths,1);
            boundary_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(displacements.m,&block_lengths(1),&displacements(1)); // TODO: collapse consecutive elements into blocks
            boundary_datatypes(s).Commit();}
        int ghost_indices_length=partition.ghost_indices(s).Size()+1;
        if(ghost_indices_length){
            ghost_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(1,&ghost_indices_length,&partition.ghost_indices(s).min_corner);
            ghost_datatypes(s).Commit();}}
}
//#####################################################################
// Function Find_Ghost_Regions
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Find_Ghost_Regions(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries)
{
    // Find which columns we need from each of the other procs
    columns_to_receive.Resize(proc_column_index_boundaries.m);
    columns_to_send.Resize(proc_column_index_boundaries.m);
    VECTOR_ND<bool> column_needed(A.n);
    int row_index=A.offsets(1);
    for(int row=1;row<=A.n;row++){ // First out which columns of A we actually have stuff in
        const int end=A.offsets(row+1);
        for(;row_index<end;row_index++){
            column_needed(A.A(row_index).j)=true;}}
    int my_rank=comm.Get_rank();ARRAY<int> temp_indices;temp_indices.Preallocate(proc_column_index_boundaries(1).y-proc_column_index_boundaries(1).x);
    for(int node_rank=1;node_rank<=proc_column_index_boundaries.m;node_rank++){
        if(node_rank-1!=my_rank){
            temp_indices.Remove_All();
            for(int column_index=proc_column_index_boundaries(node_rank).x;column_index<=proc_column_index_boundaries(node_rank).y;column_index++)
                if(column_needed(column_index)) temp_indices.Append(column_index);
            temp_indices.Compact();
            columns_to_receive(node_rank)=temp_indices;}}

    // Send this to each of the other procs
    ARRAY<ARRAY<char> > send_buffers(columns_to_send.m);ARRAY<MPI::Request> requests;
    int tag=1; // TODO change this to an actual tag value
    for(int node_rank=1;node_rank<=columns_to_receive.m;node_rank++){
        if(node_rank-1!=my_rank){
            int buffer_size=1+MPI_UTILITIES::Pack_Size(columns_to_receive(node_rank),comm);
            send_buffers(node_rank).Resize(buffer_size);int position=0;
            MPI_UTILITIES::Pack(columns_to_receive(node_rank),send_buffers(node_rank),position,comm);
            requests.Append(comm.Isend(&(send_buffers(node_rank)(1)),position,MPI::PACKED,node_rank-1,tag));}}
    // Receive a list from each of the other procs, and store this list
    for(int node_rank=1;node_rank<=columns_to_receive.m-1;node_rank++){
        MPI::Status status;
        comm.Probe(MPI::ANY_SOURCE,tag,status);
        int source=status.Get_source();
        ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
        requests.Append(comm.Irecv(&buffer(1),buffer.m,MPI::PACKED,source,tag));
        MPI_UTILITIES::Unpack(columns_to_send(source+1),buffer,position,comm);}
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Find_Ghost_Regions
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Find_Ghost_Regions_Threaded(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries)
{
    // Find which columns we need from each of the other procs
    columns_to_receive.Resize(proc_column_index_boundaries.m);
    columns_to_send.Resize(proc_column_index_boundaries.m);
    VECTOR_ND<bool> column_needed(A.n);
    int row_index=A.offsets(1);
    for(int row=1;row<=A.n;row++){ // First out which columns of A we actually have stuff in
        const int end=A.offsets(row+1);
        for(;row_index<end;row_index++){
            column_needed(A.A(row_index).j)=true;}}
    int my_rank=thread_grid->rank;ARRAY<int> temp_indices;temp_indices.Preallocate(proc_column_index_boundaries(1).y-proc_column_index_boundaries(1).x);
    for(int node_rank=1;node_rank<=proc_column_index_boundaries.m;node_rank++){
        if(node_rank-1!=my_rank){
            temp_indices.Remove_All();
            for(int column_index=proc_column_index_boundaries(node_rank).x;column_index<=proc_column_index_boundaries(node_rank).y;column_index++)
                if(column_needed(column_index)) temp_indices.Append(column_index);
            temp_indices.Compact();
            columns_to_receive(node_rank)=temp_indices;}}

    for(int node_rank=1;node_rank<=columns_to_receive.m;node_rank++){
        if(node_rank-1!=my_rank){
            int buffer_size=sizeof(int)+sizeof(T)*columns_to_receive(node_rank).m;
            THREAD_PACKAGE pack(buffer_size);int position=0;pack.send_tid=my_rank;pack.recv_tid=node_rank-1;
            *(int*)(&pack.buffer(position+1))=columns_to_receive(node_rank).m;position+=sizeof(int);
            for(int i=1;i<=columns_to_receive(node_rank).m;i++){*(int*)(&pack.buffer(position+1))=columns_to_receive(node_rank)(i);position+=sizeof(int);}
            pthread_mutex_lock(thread_grid->lock);
            thread_grid->buffers.Append(pack);
            pthread_mutex_unlock(thread_grid->lock);}}
    // Receive a list from each of the other procs, and store this list
    pthread_barrier_wait(thread_grid->barr);
    for(int buf=1;buf<=thread_grid->buffers.m;buf++){if(thread_grid->buffers(buf).recv_tid!=my_rank) continue;
        int source=thread_grid->buffers(buf).send_tid;int position=0;
        columns_to_send(source+1).Resize(*(int*)(&thread_grid->buffers(buf).buffer(position+1)));position+=sizeof(int);
        for(int i=1;i<=columns_to_send(source+1).m;i++){columns_to_send(source+1)(i)=*(int*)(&thread_grid->buffers(buf).buffer(position+1));position+=sizeof(int);}}
    pthread_barrier_wait(thread_grid->barr);
    if(thread_grid->tid==1) thread_grid->buffers.m=0;
    pthread_barrier_wait(thread_grid->barr);
}
//#####################################################################
// Function Fill_Ghost_Cells_Far
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Fill_Ghost_Cells_Far(VECTOR_ND<T>& x)
{
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;requests.Preallocate(2*columns_to_send.m);
    int my_rank=comm.Get_rank();int tag=1;
    // Send out the column values that we owe other people
    for(int node_rank=1;node_rank<=columns_to_send.m;node_rank++){
        if(node_rank-1!=my_rank){
            // First build the array of column values wanted
            ARRAY<T> send_array(columns_to_send(node_rank).m);
            for(int i=1;i<=columns_to_send(node_rank).m;i++)
                send_array(i)=x(columns_to_send(node_rank)(i));
            MPI_PACKAGE package(send_array);
            packages.Append(package);requests.Append(package.Isend(comm,node_rank-1,tag));}}
    // Receive the column values that others owe us
    ARRAY<ARRAY<T> > columns_to_receive_values(columns_to_receive.m);
    for(int node_rank=1;node_rank<=columns_to_receive.m;node_rank++){
        if(node_rank-1!=my_rank){
            // First build the array of column values we will receive
            columns_to_receive_values(node_rank).Resize(columns_to_receive(node_rank).m);
            MPI_PACKAGE package(columns_to_receive_values(node_rank));
            packages.Append(package);requests.Append(package.Irecv(comm,node_rank-1,tag));}}
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);

    // For the ones we received, stick them into x
    for(int node_rank=1;node_rank<=columns_to_receive.m;node_rank++){
        if(node_rank-1!=my_rank)
            for(int i=1;i<=columns_to_receive(node_rank).m;i++)
                x(columns_to_receive(node_rank)(i))=columns_to_receive_values(node_rank)(i);}
}
//#####################################################################
// Function Fill_Ghost_Cells_Far
//#####################################################################
template<class T_GRID> void PCG_SPARSE_MPI<T_GRID>::
Fill_Ghost_Cells_Threaded(VECTOR_ND<T>& x)
{
    int my_rank=thread_grid->rank;
    // Send out the column values that we owe other people
    for(int node_rank=1;node_rank<=columns_to_send.m;node_rank++){
        if(node_rank-1!=my_rank){
            // First build the array of column values wanted
            THREAD_PACKAGE pack(columns_to_send(node_rank).m*sizeof(T));pack.send_tid=my_rank;pack.recv_tid=node_rank-1;int position=0;
            for(int i=1;i<=columns_to_send(node_rank).m;i++){*(T*)(pack.buffer(position+1))=x(columns_to_send(node_rank)(i));position+=sizeof(T);}
            pthread_mutex_lock(thread_grid->lock);
            thread_grid->buffers.Append(pack);
            pthread_mutex_unlock(thread_grid->lock);}}
    // Receive the column values that others owe us
    ARRAY<ARRAY<T> > columns_to_receive_values(columns_to_receive.m);
    pthread_barrier_wait(thread_grid->barr);
    for(int buf=1;buf<=thread_grid->buffers.m;buf++){if(my_rank!=thread_grid->buffers(buf).recv_tid) continue;
        // First build the array of column values we will receive
        int node_rank=thread_grid->buffers(buf).send_tid;columns_to_receive_values(node_rank).Resize(columns_to_receive(node_rank).m);int position=0;
        for(int i=1;i<=columns_to_receive_values(node_rank).m;i++){columns_to_receive_values(node_rank)(i)=*(T*)(&thread_grid->buffers(buf).buffer(position+1));position+=sizeof(T);}}
    pthread_barrier_wait(thread_grid->barr);
    if(thread_grid->tid==1) thread_grid->buffers.m=0;
    pthread_barrier_wait(thread_grid->barr);

    // For the ones we received, stick them into x
    for(int node_rank=1;node_rank<=columns_to_receive.m;node_rank++){
        if(node_rank-1!=my_rank)
            for(int i=1;i<=columns_to_receive(node_rank).m;i++)
                x(columns_to_receive(node_rank)(i))=columns_to_receive_values(node_rank)(i);}
}
//#####################################################################
template class PCG_SPARSE_MPI<GRID<VECTOR<float,1> > >;
template class PCG_SPARSE_MPI<GRID<VECTOR<float,2> > >;
template class PCG_SPARSE_MPI<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PCG_SPARSE_MPI<GRID<VECTOR<double,1> > >;
template class PCG_SPARSE_MPI<GRID<VECTOR<double,2> > >;
template class PCG_SPARSE_MPI<GRID<VECTOR<double,3> > >;
#endif
#endif
