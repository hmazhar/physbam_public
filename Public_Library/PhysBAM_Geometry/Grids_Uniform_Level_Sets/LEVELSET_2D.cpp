//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LEVELSET_2D<T_GRID>::
LEVELSET_2D(T_GRID& grid_input,ARRAY<T,VECTOR<int,2> >& phi_input,const int number_of_ghost_cells_input)
    :LEVELSET_UNIFORM<T_GRID>(grid_input,phi_input,number_of_ghost_cells_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LEVELSET_2D<T_GRID>::
~LEVELSET_2D()
{}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T_GRID> SYMMETRIC_MATRIX<typename T_GRID::SCALAR,2> LEVELSET_2D<T_GRID>::
Hessian(const VECTOR<T,2>& X) const
{
    T one_over_dx=1/grid.dX.x,one_over_dy=1/grid.dX.y;T two_phi_center=2*Phi(X);
    T phi_xx=(Phi(VECTOR<T,2>(X.x+grid.dX.x,X.y))-two_phi_center+Phi(VECTOR<T,2>(X.x-grid.dX.x,X.y)))*sqr(one_over_dx),
       phi_yy=(Phi(VECTOR<T,2>(X.x,X.y+grid.dX.y))-two_phi_center+Phi(VECTOR<T,2>(X.x,X.y-grid.dX.y)))*sqr(one_over_dy),
       phi_xy=(Phi(VECTOR<T,2>(X.x+grid.dX.x,X.y+grid.dX.y))-Phi(VECTOR<T,2>(X.x+grid.dX.x,X.y-grid.dX.y))
                  -Phi(VECTOR<T,2>(X.x-grid.dX.x,X.y+grid.dX.y))+Phi(VECTOR<T,2>(X.x-grid.dX.x,X.y-grid.dX.y)))*(T).25*one_over_dx*one_over_dy;
    return SYMMETRIC_MATRIX<T,2>(phi_xx,phi_xy,phi_yy);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,1> LEVELSET_2D<T_GRID>::
Principal_Curvatures(const VECTOR<T,2>& X) const
{
    VECTOR<T,2> grad_phi=VECTOR<T,2>((Phi(VECTOR<T,2>(X.x+grid.dX.x,X.y))-Phi(VECTOR<T,2>(X.x-grid.dX.x,X.y)))/(2*grid.dX.x),
                                     (Phi(VECTOR<T,2>(X.x,X.y+grid.dX.y))-Phi(VECTOR<T,2>(X.x,X.y-grid.dX.y)))/(2*grid.dX.y));
    VECTOR<T,2> tangent=grad_phi.Perpendicular();T grad_phi_magnitude=tangent.Normalize();
    T curvature=VECTOR<T,2>::Dot_Product(tangent,Hessian(X)*tangent)/grad_phi_magnitude;
    return VECTOR<T,1>(curvature);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
// note that sqrt(phix^2+phiy^2)=1 if it's a distance function
template<class T_GRID> void LEVELSET_2D<T_GRID>::
Compute_Normals(const T time)
{
    T one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y);
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
        
    if(!normals) normals=new ARRAY<VECTOR<T,2> ,VECTOR<int,2> >(grid.Domain_Indices(ghost_cells-1));
    for(int i=normals->domain.min_corner.x;i<=normals->domain.max_corner.x;i++) for(int j=normals->domain.min_corner.y;j<=normals->domain.max_corner.y;j++){
        (*normals)(i,j)=VECTOR<T,2>((phi_ghost(i+1,j)-phi_ghost(i-1,j))*one_over_two_dx,(phi_ghost(i,j+1)-phi_ghost(i,j-1))*one_over_two_dy).Normalized();}
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
// kappa = - DIV(normal), negative for negative phi inside, positive for positive phi inside, sqrt(phix^2+phiy^2)=1 for a distance function
template<class T_GRID> void LEVELSET_2D<T_GRID>::
Compute_Curvature(const T time)
{
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);

    if(!curvature) curvature=new ARRAY<T,VECTOR<int,2> >(grid.Domain_Indices(2));
    for(CELL_ITERATOR iterator(grid,2);iterator.Valid();iterator.Next())
        (*curvature)(iterator.Cell_Index())=Compute_Curvature(phi_ghost,iterator.Cell_Index());
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_2D<T_GRID>::
Compute_Curvature(const ARRAY<T,VECTOR<int,2> >& phi_input,const VECTOR<int,2>& index) const
{
    T one_over_two_dx=(T).5*grid.one_over_dX.x,one_over_two_dy=(T).5*grid.one_over_dX.y;
    T one_over_dx_squared=sqr(grid.one_over_dX.x),one_over_dy_squared=sqr(grid.one_over_dX.y),one_over_four_dx_dy=(T).25*grid.one_over_dX.x*grid.one_over_dX.y;
    T max_curvature=1/grid.min_dX; // max resolution

    int i=index.x,j=index.y;
    T phix=(phi_input(i+1,j)-phi_input(i-1,j))*one_over_two_dx,phixx=(phi_input(i+1,j)-2*phi_input(i,j)+phi_input(i-1,j))*one_over_dx_squared,
        phiy=(phi_input(i,j+1)-phi_input(i,j-1))*one_over_two_dy,phiyy=(phi_input(i,j+1)-2*phi_input(i,j)+phi_input(i,j-1))*one_over_dy_squared,
        phixy=(phi_input(i+1,j+1)-phi_input(i+1,j-1)-phi_input(i-1,j+1)+phi_input(i-1,j-1))*one_over_four_dx_dy;
    T denominator=sqrt(sqr(phix)+sqr(phiy)),curvature;
    if(denominator >= small_number) curvature=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx)/cube(denominator);else curvature=LEVELSET_UTILITIES<T>::Sign(phi_input(i,j))*max_curvature;
    return minmag(curvature,LEVELSET_UTILITIES<T>::Sign(curvature)*max_curvature);
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_2D<T_GRID>::
Compute_Curvature(const VECTOR<T,2>& location) const
{
    TV l2=(location-(T).5*grid.dX-grid.domain.min_corner)*grid.one_over_dX+1;
    TV_INT cell(floor(l2));

    TV w(l2-TV(cell));
    T k00=Compute_Curvature(phi,cell);
    cell.y++;
    T k01=Compute_Curvature(phi,cell);
    cell.x++;
    T k11=Compute_Curvature(phi,cell);
    cell.y--;
    T k10=Compute_Curvature(phi,cell);
    return LINEAR_INTERPOLATION<T,T>::Bilinear(k00,k10,k01,k11,w);
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class T_GRID> void LEVELSET_2D<T_GRID>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    if(!recompute_if_exists && cell_range) return;
    if(!cell_range) cell_range=new ARRAY<RANGE<VECTOR<T,1> >,VECTOR<int,2> >(phi.domain.min_corner.x,phi.domain.max_corner.x-1,phi.domain.min_corner.y,phi.domain.max_corner.y-1);
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x-1;i++) for(int j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y-1;j++){
        T phi1=phi(i,j),phi2=phi(i,j+1),phi3=phi(i+1,j),phi4=phi(i+1,j+1);
        (*cell_range)(i,j)=RANGE<VECTOR<T,1> >(min(phi1,phi2,phi3,phi4),max(phi1,phi2,phi3,phi4));}
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void LEVELSET_2D<T_GRID>::
Fast_Marching_Method(const T time,const T stopping_distance,const ARRAY<VECTOR<int,2> >* seed_indices,const bool add_seed_indices_for_ghost_cells)
{       
    Get_Signed_Distance_Using_FMM(phi,time,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class T_GRID> void LEVELSET_2D<T_GRID>::
Get_Signed_Distance_Using_FMM(ARRAY<T,VECTOR<int,2> >& signed_distance,const T time,const T stopping_distance,const ARRAY<VECTOR<int,2> >* seed_indices,const bool add_seed_indices_for_ghost_cells)
{       
    const int ghost_cells=2*number_of_ghost_cells+1;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<T_GRID> fmm(*this,ghost_cells);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
    ARRAY<T,VECTOR<int,2> >::Get(signed_distance,phi_ghost);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
}
//#####################################################################
// Function Fast_Marching_Method_Outside_Band
//#####################################################################
template<class T_GRID> void LEVELSET_2D<T_GRID>::
Fast_Marching_Method_Outside_Band(const T half_band_width,const T time,const T stopping_distance)
{              
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<T_GRID> fmm(*this,ghost_cells);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance);
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) if(abs(phi_ghost(i,j)) > half_band_width) phi(i,j)=phi_ghost(i,j);
    boundary->Apply_Boundary_Condition(grid,phi,time);
}
//#####################################################################
// Function Approximate_Length
//#####################################################################
// calculates the approximate perimeter using delta functions
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_2D<T_GRID>::
Approximate_Length(const T interface_thickness,const T time) const
{
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    T interface_half_width=interface_thickness*grid.dX.Max()/2,one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y),length=0;
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++)
        length+=(LEVELSET_UTILITIES<T>::Delta(phi_ghost(i,j),interface_half_width)*sqrt(sqr((phi_ghost(i+1,j)-phi_ghost(i-1,j))*one_over_two_dx)+sqr((phi_ghost(i,j+1)-phi_ghost(i,j-1))*one_over_two_dy)));
    return length*grid.dX.x*grid.dX.y;
}
//#####################################################################
template class LEVELSET_2D<GRID<VECTOR<float,2> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_2D<GRID<VECTOR<double,2> > >;
#endif
