//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_1D_HELPER.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> LEVELSET_1D<T>::
LEVELSET_1D(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,1> >& phi_input,const int number_of_ghost_cells_input)
    :LEVELSET_UNIFORM<GRID<TV> >(grid_input,phi_input,number_of_ghost_cells_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> LEVELSET_1D<T>::
~LEVELSET_1D()
{}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> MATRIX<T,1> LEVELSET_1D<T>::
Hessian(const VECTOR<T,1>& X) const
{
    return MATRIX<T,1>((Phi(X+grid.dX)-2*Phi(X)+Phi(X-grid.dX))/sqr(grid.dX.x));
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,0> LEVELSET_1D<T>::
Principal_Curvatures(const VECTOR<T,1>& X) const
{
    return VECTOR<T,0>(); // not much curvature in 1D
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
// note that abs(phix)=1 if it's a distance function
template<class T> void LEVELSET_1D<T>::
Compute_Normals(const T time)
{
    T one_over_two_dx=1/(2*grid.dX.x);
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
        
    if(!normals) normals=new ARRAY<VECTOR<T,1> ,VECTOR<int,1> >(grid.Domain_Indices(ghost_cells-1));
    for(int i=normals->domain.min_corner.x;i<=normals->domain.max_corner.x;i++)(*normals)(i)=VECTOR<T,1>((phi_ghost(i+1)-phi_ghost(i-1))*one_over_two_dx).Normalized();
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T> void LEVELSET_1D<T>::
Compute_Curvature(const T time)
{      
    if(!curvature) curvature=new ARRAY<T,VECTOR<int,1> >(grid.counts.x,2);
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class T> void LEVELSET_1D<T>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    if(!recompute_if_exists && cell_range) return;
    if(!cell_range) cell_range=new ARRAY<RANGE<VECTOR<T,1> >,VECTOR<int,1> >(phi.domain.min_corner.x,phi.domain.max_corner.x-1);
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x-1;i++){
        T phi1=phi(i),phi2=phi(i+1);
        (*cell_range)(i)=RANGE<VECTOR<T,1> >(min(phi1,phi2),max(phi1,phi2));}
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T> void LEVELSET_1D<T>::
Fast_Marching_Method(const T time,const T stopping_distance,const ARRAY<VECTOR<int,1> >* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    Get_Signed_Distance_Using_FMM(phi,time,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class T> void LEVELSET_1D<T>::
Get_Signed_Distance_Using_FMM(ARRAY<T,VECTOR<int,1> >& signed_distance,const T time,const T stopping_distance,const ARRAY<VECTOR<int,1> >* seed_indices,
        const bool add_seed_indices_for_ghost_cells)
{       
    const int ghost_cells=2*number_of_ghost_cells+1;
    ARRAY<T,VECTOR<int,1> > phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(*this,ghost_cells);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
    ARRAY<T,VECTOR<int,1> >::Get(signed_distance,phi_ghost);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
}
//#####################################################################
// Function Fast_Marching_Method_Outside_Band
//#####################################################################
template<class T> void LEVELSET_1D<T>::
Fast_Marching_Method_Outside_Band(const T half_band_width,const T time,const T stopping_distance)
{              
    int m=grid.counts.x;
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,VECTOR<int,1> > phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(*this,ghost_cells);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance);
    for(int i=1;i<=m;i++) if(abs(phi_ghost(i)) > half_band_width) phi(i)=phi_ghost(i);
    boundary->Apply_Boundary_Condition(grid,phi,time);
}
//#####################################################################
template class LEVELSET_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_1D<double>;
#endif
