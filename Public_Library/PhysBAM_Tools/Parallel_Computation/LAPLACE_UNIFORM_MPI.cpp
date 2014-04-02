//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Frank Losasso, Michael Lentine, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_UNIFORM_MPI.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Data_Structures/GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#endif
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,1>&)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,1> >(1,m),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,1> >(0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,1> >(m+1,m+1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,1> >(1,1),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,1> >(m,m),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x,n=local_grid.counts.y;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,2> >(1,m,1,n),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,2> >(0,0,1,n),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,2> >(m+1,m+1,1,n),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(3,RANGE<VECTOR<int,2> >(1,m,0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(4,RANGE<VECTOR<int,2> >(1,m,n+1,n+1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,2> >(1,1,1,n),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,2> >(m,m,1,n),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(3,RANGE<VECTOR<int,2> >(1,m,1,1),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(4,RANGE<VECTOR<int,2> >(1,m,n,n),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x,n=local_grid.counts.y,mn=local_grid.counts.z;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,3> >(1,m,1,n,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,3> >(0,0,1,n,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,3> >(m+1,m+1,1,n,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(3,RANGE<VECTOR<int,3> >(1,m,0,0,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(4,RANGE<VECTOR<int,3> >(1,m,n+1,n+1,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(5,RANGE<VECTOR<int,3> >(1,m,1,n,0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(6,RANGE<VECTOR<int,3> >(1,m,1,n,mn+1,mn+1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,3> >(1,1,1,n,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,3> >(m,m,1,n,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(3,RANGE<VECTOR<int,3> >(1,m,1,1,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(4,RANGE<VECTOR<int,3> >(1,m,n,n,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(5,RANGE<VECTOR<int,3> >(1,m,1,n,1,1),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(6,RANGE<VECTOR<int,3> >(1,m,1,n,mn,mn),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::
Find_Matrix_Indices_Threaded(ARRAY<RANGE<TV_INT> >& domains,ARRAY<ARRAY<INTERVAL<int> > >& interior_indices,ARRAY<ARRAY<ARRAY<INTERVAL<int> > > >& ghost_indices,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,LAPLACE_UNIFORM<T_GRID>* laplace)
{
    assert(local_grid.Is_MAC_Grid());
    //interior mpi cells, interior thread cells
    for(int color=1;color<=filled_region_ranks.m;color++) partitions(color).interior_indices.min_corner=filled_region_cell_count(color)+1;
    for(int i=1;i<=domains.m;i++){
        RANGE<TV_INT> interior_domain(domains(i));
        interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
        for(int color=1;color<=interior_indices.m;color++) interior_indices(color)(i).min_corner=filled_region_cell_count(color)+1;
        laplace->Compute_Matrix_Indices(interior_domain,filled_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
        for(int color=1;color<=interior_indices.m;color++) interior_indices(color)(i).max_corner=filled_region_cell_count(color);}
    for(int color=1;color<=filled_region_ranks.m;color++) partitions(color).interior_indices.max_corner=filled_region_cell_count(color);
    //boundary mpi cells
    for(int axis=1;axis<=TV::dimension;axis++) for(int side=1;side<=2;side++){int s=(axis-1)*2+side;
        for(int color=1;color<=filled_region_ranks.m;color++) partitions(color).ghost_indices(s).min_corner=filled_region_cell_count(color)+1;
        RANGE<TV_INT> exterior_domain(local_grid.Domain_Indices(1));
        for(int axis2=axis+1;axis2<=TV::dimension;axis2++){exterior_domain.min_corner(axis2)++;exterior_domain.max_corner(axis2)--;}
        if(side==1) exterior_domain.max_corner(axis)=exterior_domain.min_corner(axis);
        else exterior_domain.min_corner(axis)=exterior_domain.max_corner(axis);
        for(int i=1;i<=domains.m;i++){
            RANGE<TV_INT> interior_domain(domains(i));
            interior_domain.max_corner-=TV_INT::All_Ones_Vector();for(int axis=1;axis<=TV_INT::dimension;axis++) if(interior_domain.max_corner(axis)==local_grid.Domain_Indices().max_corner(axis)) interior_domain.max_corner(axis)++;
            interior_domain.min_corner+=TV_INT::All_Ones_Vector();for(int axis=1;axis<=TV_INT::dimension;axis++) if(interior_domain.min_corner(axis)==local_grid.Domain_Indices().min_corner(axis)) interior_domain.min_corner(axis)--;
            for(int color=1;color<=interior_indices.m;color++) ghost_indices(color)(i)(s).min_corner=filled_region_cell_count(color)+1;
            laplace->Compute_Matrix_Indices(RANGE<TV_INT>::Intersect(exterior_domain,interior_domain),filled_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
            for(int color=1;color<=interior_indices.m;color++) ghost_indices(color)(i)(s).max_corner=filled_region_cell_count(color);}
        for(int color=1;color<=filled_region_ranks.m;color++) partitions(color).ghost_indices(s).max_corner=filled_region_cell_count(color);
        Find_Boundary_Indices_In_Region(s,exterior_domain,cell_index_to_matrix_index);}
}
//#####################################################################
// Function Find_Matrix_Indices_In_Region
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::
Find_Matrix_Indices_In_Region(const int region_index,const RANGE<TV_INT>& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,
    ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array)
{
    if(region_index) for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).min_corner=filled_region_cell_count(color)+1;
    else for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).interior_indices.min_corner=filled_region_cell_count(color)+1;
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){TV_INT c=iterator.Cell_Index();
        int color=filled_region_colors(c);if(color<1 || (!filled_region_touches_dirichlet(color)&&!solve_neumann_regions)) continue;
        int new_index=++filled_region_cell_count(color);cell_index_to_matrix_index(c)=new_index;
        matrix_index_to_cell_index_array(color)(new_index)=c;}
    if(region_index) for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).max_corner=filled_region_cell_count(color);
    else for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).interior_indices.max_corner=filled_region_cell_count(color);
}
//#####################################################################
// Function Find_Boundary_Indices_In_Region
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::
Find_Boundary_Indices_In_Region(const int side,const RANGE<TV_INT>& region,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    int axis=(side-1)/2+1,cell_side=side&1;
    RANGE<TV_INT> face_region=region;if(!cell_side) face_region+=TV_INT::Axis_Vector(axis);
    // count boundary indices
    ARRAY<int> counts(partitions.m);
    TV_INT face_offset=cell_side?TV_INT():TV_INT::Axis_Vector(axis);
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next())if(!psi_N.Component(axis)(iterator.Cell_Index()+face_offset)){
        int color=filled_region_colors(iterator.Cell_Index());
        if(counts.Valid_Index(color)) counts(color)++;}
    // fill boundary indices
    for(int color=1;color<=partitions.m;color++)partitions(color).boundary_indices(side).Resize(counts(color));
    ARRAYS_COMPUTATIONS::Fill(counts,0);
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next())if(!psi_N.Component(axis)(iterator.Cell_Index()+face_offset)){
        int color=filled_region_colors(iterator.Cell_Index());
        if(counts.Valid_Index(color)) partitions(color).boundary_indices(side)(++counts(color))=cell_index_to_matrix_index(iterator.Cell_Index());}
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const VECTOR<T,1>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const VECTOR<T,2>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const VECTOR<T,3>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> void LAPLACE_UNIFORM_MPI<T_GRID>::Find_Matrix_Indices_Threaded(ARRAY<RANGE<TV_INT> >&,ARRAY<ARRAY<INTERVAL<int> > >&,ARRAY<ARRAY<ARRAY<INTERVAL<int> > > >&,ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,LAPLACE_UNIFORM<T_GRID>*){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################

#endif

//#####################################################################
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER(T_GRID) \
    template void LAPLACE_UNIFORM_MPI<T_GRID >::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const TV&); \
    template void LAPLACE_UNIFORM_MPI<T_GRID >::Find_Matrix_Indices_Threaded(ARRAY<RANGE<TV_INT> >&,ARRAY<ARRAY<INTERVAL<int> > >&,ARRAY<ARRAY<ARRAY<INTERVAL<int> > > >&,ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,LAPLACE_UNIFORM<T_GRID>*);
template LAPLACE_UNIFORM_MPI<GRID<VECTOR<float,1> > >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<float,1> >) >&);
template LAPLACE_UNIFORM_MPI<GRID<VECTOR<float,2> > >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<float,2> >) >&);
template LAPLACE_UNIFORM_MPI<GRID<VECTOR<float,3> > >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<float,3> >) >&);
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template LAPLACE_UNIFORM_MPI<GRID<VECTOR<double,1> > >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<double,1> >) >&);
template LAPLACE_UNIFORM_MPI<GRID<VECTOR<double,2> > >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<double,2> >) >&);
template LAPLACE_UNIFORM_MPI<GRID<VECTOR<double,3> > >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<double,3> >) >&);
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >));
#endif
