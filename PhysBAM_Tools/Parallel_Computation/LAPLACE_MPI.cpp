//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_PDE_Linear/LAPLACE_RLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#endif
#include <PhysBAM_Tools/Parallel_Computation/FLOOD_FILL_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI_THREADED.h>
#endif
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LAPLACE_MPI<T_GRID>::
LAPLACE_MPI(T_LAPLACE& laplace)
    :mpi_grid(laplace.mpi_grid),local_grid(laplace.grid),local_pcg(laplace.pcg),number_of_regions(laplace.number_of_regions),number_of_global_regions(-1),filled_region_colors(laplace.filled_region_colors),
    filled_region_touches_dirichlet(laplace.filled_region_touches_dirichlet),solve_neumann_regions(laplace.solve_neumann_regions),psi_N(laplace.psi_N)
{
    groups=new ARRAY<MPI::Group>;
    communicators=new ARRAY<MPI::Intracomm>;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LAPLACE_MPI<T_GRID>::
~LAPLACE_MPI()
{
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(*groups);delete groups;
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(*communicators);delete communicators;
}
//#####################################################################
// Function Synchronize_Solution_Regions
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Synchronize_Solution_Regions()
{
    FLOOD_FILL_MPI<T_GRID> flood_fill(*mpi_grid,local_grid,psi_N,number_of_regions,filled_region_colors,filled_region_ranks,&filled_region_touches_dirichlet);
    number_of_global_regions=flood_fill.Synchronize_Colors();

    // allocate communicators for each color
    ARRAY<MPI::Group> new_groups(filled_region_ranks.m);
    ARRAY<MPI::Intracomm> new_communicators(filled_region_ranks.m);
    for(int color=1;color<=filled_region_ranks.m;color++){
        new_groups(color)=mpi_grid->group->Incl(filled_region_ranks(color).m,&filled_region_ranks(color)(1));
        int i;
        for(i=1;i<=groups->m;i++)if((*groups)(i)!=MPI::GROUP_NULL && MPI::Group::Compare((*groups)(i),new_groups(color))==MPI::IDENT){
           new_communicators(color)=(*communicators)(i);(*communicators)(i)=MPI::COMM_NULL;(*groups)(i).Free();break;}
        if(i>groups->m) new_communicators(color)=mpi_grid->comm->Create(new_groups(color));} // NOTE: returns null communicator for processes not in the group!
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(*groups);MPI_UTILITIES::Free_Elements_And_Clean_Memory(*communicators);
    groups->Exchange(new_groups);communicators->Exchange(new_communicators);

    // allocate partitions and compute neighbor ranks
    partitions.Resize(filled_region_ranks.m);
    for(int color=1;color<=filled_region_ranks.m;color++){
        partitions(color).Set_Number_Of_Sides(T_PARALLEL_GRID::number_of_faces_per_cell);
        for(int s=1;s<=T_PARALLEL_GRID::number_of_faces_per_cell;s++){
            int global_rank=mpi_grid->side_neighbor_ranks(s);
            if(global_rank==MPI::PROC_NULL) partitions(color).neighbor_ranks(s)=MPI::PROC_NULL;
            else MPI::Group::Translate_ranks(*mpi_grid->group,1,&global_rank,(*groups)(color),&partitions(color).neighbor_ranks(s));}}
}
//#####################################################################
// Function Update_Solution_Regions_For_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Update_Solution_Regions_For_Solid_Fluid_Coupling(const T_MPI_GRID& mpi_grid)
{
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(*groups);MPI_UTILITIES::Free_Elements_And_Clean_Memory(*communicators);
    groups->Resize(1);communicators->Resize(1);
    (*communicators)(1)=mpi_grid.comm->Dup();(*groups)(1)=*new MPI::Group((*communicators)(1).Get_group());
}
//#####################################################################
// Function Update_Solution_Regions_For_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> int LAPLACE_MPI<T_GRID>::
Get_Total_Number_Of_Threads(const int input,const int color)
{
    int output;
    MPI_UTILITIES::Reduce(input,output,MPI::SUM,(*communicators)(color));
    return output;
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Solve_Threaded(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,ARRAY<INTERVAL<int> >& interior_indices,ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance,const int color,const int multi_proc_mode)
{
    if(color>filled_region_ranks.m){
        if(multi_proc_mode){local_pcg_threaded->Solve(domain,domain_index,interior_indices,ghost_indices,A,x,b,tolerance);return;}
        else{local_pcg.Solve(A,x,b,q,s,r,k,z,tolerance);return;}}
    else{
        PCG_SPARSE_MPI_THREADED<TV> pcg_mpi(*local_pcg_threaded,(*communicators)(color),partitions(color));
        assert(Use_Parallel_Solve());
        pcg_mpi.Solve(domain,domain_index,interior_indices,ghost_indices,A,x,b,tolerance);}
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance,const int color)
{
    if(color>filled_region_ranks.m){local_pcg.Solve(A,x,b,q,s,r,k,z,tolerance);return;}
    else{
        PCG_SPARSE_MPI<T_GRID> pcg_mpi(local_pcg,(*communicators)(color),partitions(color));
        pcg_mpi.thread_grid=mpi_grid->threaded_grid;
        if(Use_Parallel_Solve()) pcg_mpi.Parallel_Solve(A,x,b,tolerance);
        else pcg_mpi.Serial_Solve(A,x,b,q,s,r,k,z,1234,tolerance);}
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_MPI<T_GRID>::
Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance,const int color,const ARRAY<VECTOR<int,2> >& global_column_index_boundaries)
{
    SPARSE_MATRIX_PARTITION temp_partition;
    PCG_SPARSE_MPI<T_GRID> pcg_mpi(local_pcg,(*communicators)(color),temp_partition);
    pcg_mpi.thread_grid=mpi_grid->threaded_grid;
    pcg_mpi.Parallel_Solve(A,x,b,global_column_index_boundaries,tolerance,true);
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID> LAPLACE_MPI<T_GRID>::LAPLACE_MPI(T_LAPLACE& laplace):mpi_grid(laplace.mpi_grid),local_grid(laplace.grid),local_pcg(laplace.pcg),
    number_of_regions(laplace.number_of_regions),filled_region_colors(laplace.filled_region_colors),filled_region_touches_dirichlet(laplace.filled_region_touches_dirichlet),
    solve_neumann_regions(laplace.solve_neumann_regions),psi_N(laplace.psi_N),groups(0),communicators(0){}
template<class T_GRID> LAPLACE_MPI<T_GRID>::~LAPLACE_MPI(){}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Synchronize_Solution_Regions(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Update_Solution_Regions_For_Solid_Fluid_Coupling(const T_MPI_GRID& mpi_grid){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Solve(SPARSE_MATRIX_FLAT_NXN<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,
    const T,const int){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance,const int color,
    const ARRAY<VECTOR<int,2> >& global_column_index_boundaries){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_MPI<T_GRID>::Solve_Threaded(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,ARRAY<INTERVAL<int> >& interior_indices,
    ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,
    VECTOR_ND<T>& z,const T tolerance,const int color,const int multi_proc_mode)
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################

#endif

template<class T_GRID> bool& LAPLACE_MPI<T_GRID>::
Use_Parallel_Solve()
{
    static bool use_parallel_solve=true;
    return use_parallel_solve;
}

//#####################################################################
template class LAPLACE_MPI<GRID<VECTOR<float,1> > >;
template class LAPLACE_MPI<GRID<VECTOR<float,2> > >;
template class LAPLACE_MPI<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_MPI<GRID<VECTOR<double,1> > >;
template class LAPLACE_MPI<GRID<VECTOR<double,2> > >;
template class LAPLACE_MPI<GRID<VECTOR<double,3> > >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class LAPLACE_MPI<RLE_GRID_2D<float> >;
template class LAPLACE_MPI<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_MPI<RLE_GRID_2D<double> >;
template class LAPLACE_MPI<RLE_GRID_3D<double> >;
#endif
#endif
