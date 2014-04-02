//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CHILDREN.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> LEVELSET_OCTREE<T>::
Hessian(const TV& X,const ARRAY<T>* phi_nodes) const
{
    TV DX(grid.Leaf_Cell(X)->DX());
    T one_over_dx=1/DX.x,one_over_dy=1/DX.y,one_over_dz=1/DX.z;T two_phi_center=2*Phi(X,phi_nodes);
    T phi_xx=(Phi(TV(X.x+DX.x,X.y,X.z),phi_nodes)-two_phi_center+Phi(TV(X.x-DX.x,X.y,X.z),phi_nodes))*sqr(one_over_dx),
        phi_yy=(Phi(TV(X.x,X.y+DX.y,X.z),phi_nodes)-two_phi_center+Phi(TV(X.x,X.y-DX.y,X.z),phi_nodes))*sqr(one_over_dy),
        phi_zz=(Phi(TV(X.x,X.y,X.z+DX.z),phi_nodes)-two_phi_center+Phi(TV(X.x,X.y,X.z-DX.z),phi_nodes))*sqr(one_over_dz),
        phi_xy=(Phi(TV(X.x+DX.x,X.y+DX.y,X.z),phi_nodes)-Phi(TV(X.x+DX.x,X.y-DX.y,X.z),phi_nodes)-Phi(TV(X.x-DX.x,X.y+DX.y,X.z),phi_nodes)+Phi(TV(X.x-DX.x,X.y-DX.y,X.z),phi_nodes))
                    *(T).25*one_over_dx*one_over_dy,
        phi_xz=(Phi(TV(X.x+DX.x,X.y,X.z+DX.z),phi_nodes)-Phi(TV(X.x+DX.x,X.y,X.z-DX.z),phi_nodes)-Phi(TV(X.x-DX.x,X.y,X.z+DX.z),phi_nodes)+Phi(TV(X.x-DX.x,X.y,X.z-DX.z),phi_nodes))
                    *(T).25*one_over_dx*one_over_dz,
        phi_yz=(Phi(TV(X.x,X.y+DX.y,X.z+DX.z),phi_nodes)-Phi(TV(X.x,X.y+DX.y,X.z-DX.z),phi_nodes)-Phi(TV(X.x,X.y-DX.y,X.z+DX.z),phi_nodes)+Phi(TV(X.x,X.y-DX.y,X.z-DX.z),phi_nodes))
                    *(T).25*one_over_dy*one_over_dz;
    return SYMMETRIC_MATRIX<T,3>(phi_xx,phi_xy,phi_xz,phi_yy,phi_yz,phi_zz);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> LEVELSET_OCTREE<T>::
Principal_Curvatures(const TV& X,const ARRAY<T>* phi_nodes) const
{
    TV DX(grid.Leaf_Cell(X)->DX());
    TV grad_phi=TV((Phi(TV(X.x+DX.x,X.y,X.z),phi_nodes)-Phi(TV(X.x-DX.x,X.y,X.z),phi_nodes))/((T)2*DX.x),
                   (Phi(TV(X.x,X.y+DX.y,X.z),phi_nodes)-Phi(TV(X.x,X.y-DX.y,X.z),phi_nodes))/((T)2*DX.y),
                   (Phi(TV(X.x,X.y,X.z+DX.z),phi_nodes)-Phi(TV(X.x,X.y,X.z-DX.z),phi_nodes))/((T)2*DX.z));
    T grad_phi_magnitude=grad_phi.Magnitude();TV N=grad_phi/grad_phi_magnitude;
    SYMMETRIC_MATRIX<T,3> P=(T)1-SYMMETRIC_MATRIX<T,3>::Outer_Product(N),M=SYMMETRIC_MATRIX<T,3>::Conjugate(P,Hessian(X,phi_nodes))/grad_phi_magnitude;
    T trace=M.Trace();
    QUADRATIC<T> quadratic(-1,trace,sqr(M(2,1))-M(1,1)*M(2,2)+sqr(M(3,1))-M(1,1)*M(3,3)+sqr(M(3,2))-M(2,2)*M(3,3));quadratic.Compute_Roots();
    if(quadratic.roots == 0) return (T).5*VECTOR<T,2>(trace,trace);
    else if(quadratic.roots == 1) return VECTOR<T,2>(quadratic.root1,quadratic.root1);
    else return VECTOR<T,2>(quadratic.root1,quadratic.root2);
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T> void LEVELSET_OCTREE<T>::
Compute_Curvature(const T time)
{
    ARRAY<T> phi_ghost(grid.number_of_cells);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);
    ARRAY<T> phi_nodes_ghost(grid.number_of_nodes);LINEAR_INTERPOLATION_DYADIC_HELPER<OCTREE_GRID<T> >::Interpolate_From_Cells_To_Nodes(grid,phi_ghost,phi_nodes_ghost);
    if(!curvature) curvature=new ARRAY<T>(grid.number_of_cells);else{curvature->Resize(grid.number_of_cells);ARRAYS_COMPUTATIONS::Fill(*curvature,0);}

    for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(grid,2);iterator.Valid();iterator.Next()){
        TV location=iterator.Location();
        TV DX=iterator.DX();

        T one_over_two_dx=1/(2*DX.x),one_over_two_dy=1/(2*DX.y),one_over_two_dz=1/(2*DX.z),one_over_dx_squared=1/sqr(DX.x),one_over_dy_squared=1/sqr(DX.y),one_over_dz_squared=1/sqr(DX.z),
            one_over_four_dx_dy=1/(4*DX.x*DX.y),one_over_four_dx_dz=1/(4*DX.x*DX.z),one_over_four_dy_dz=1/(4*DX.y*DX.z);
        T max_curvature=1/grid.Minimum_Edge_Length(); // max resolution
        const OCTREE_CELL<T>* cell=grid.Base_Cell(location);if(!cell)continue;
        T phi_curr=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y,location.z)),
            phi_mx=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y,location.z)),
            phi_px=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y,location.z)),
            phi_my=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y-DX.y,location.z)),
            phi_py=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y+DX.y,location.z)),
            phi_mz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y,location.z-DX.z)),
            phi_pz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y,location.z+DX.z)),
            phi_px_py=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y+DX.y,location.z)),
            phi_px_my=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y-DX.y,location.z)),
            phi_mx_py=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y+DX.y,location.z)),
            phi_mx_my=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y-DX.y,location.z)),
            phi_px_pz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y,location.z+DX.z)),
            phi_px_mz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y,location.z-DX.z)),
            phi_mx_pz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y,location.z+DX.z)),
            phi_mx_mz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y,location.z-DX.z)),
            phi_py_pz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y+DX.y,location.z+DX.z)),
            phi_py_mz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y+DX.y,location.z-DX.z)),
            phi_my_pz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y-DX.y,location.z+DX.z)),
            phi_my_mz=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y-DX.y,location.z-DX.z));
        T phix=(phi_px-phi_mx)*one_over_two_dx,phiy=(phi_py-phi_my)*one_over_two_dy,phiz=(phi_pz-phi_mz)*one_over_two_dz,
            phixx=(phi_px-2*phi_curr+phi_mx)*one_over_dx_squared,phiyy=(phi_py-2*phi_curr+phi_my)*one_over_dy_squared,phizz=(phi_pz-2*phi_curr+phi_mz)*one_over_dz_squared,
            phixy=(phi_px_py-phi_px_my-phi_mx_py+phi_mx_my)*one_over_four_dx_dy,phixz=(phi_px_pz-phi_px_mz-phi_mx_pz+phi_mx_mz)*one_over_four_dx_dz,
            phiyz=(phi_py_pz-phi_py_mz-phi_my_pz+phi_my_mz)*one_over_four_dy_dz;
        T denominator=sqrt(sqr(phix)+sqr(phiy)+sqr(phiz));
        T curvature_value;
        if(denominator >= small_number) 
            curvature_value=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx+
                              sqr(phix)*phizz-2*phix*phiz*phixz+sqr(phiz)*phixx+
                              sqr(phiy)*phizz-2*phiy*phiz*phiyz+sqr(phiz)*phiyy)/cube(denominator);
        else curvature_value=LEVELSET_UTILITIES<T>::Sign(phi_curr)*max_curvature;
        (*curvature)(iterator.Cell_Index())=minmag(curvature_value,LEVELSET_UTILITIES<T>::Sign(curvature_value)*max_curvature);}
}
//#####################################################################
// Function Curvature
//#####################################################################
template<class T> T LEVELSET_OCTREE<T>::
Curvature(const TV& location,const bool use_precomputed_curvatures_if_possible,const ARRAY<T>* phi_nodes) const
{
    if(curvature && use_precomputed_curvatures_if_possible) return LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,*curvature,location);
    else{ // warning: the following code is incredibly expensive
        TV DX(grid.Leaf_Cell(location)->DX());
        T one_over_two_dx=1/(2*DX.x),one_over_two_dy=1/(2*DX.y),one_over_two_dz=1/(2*DX.z),one_over_dx_squared=1/sqr(DX.x),one_over_dy_squared=1/sqr(DX.y),one_over_dz_squared=1/sqr(DX.z),
            one_over_four_dx_dy=1/(4*DX.x*DX.y),one_over_four_dx_dz=1/(4*DX.x*DX.z),one_over_four_dy_dz=1/(4*DX.y*DX.z);
        T max_curvature=1/grid.Minimum_Edge_Length(); // max resolution

        T phi_curr=Clamped_Phi(TV(location.x,location.y,location.z),phi_nodes),
            phi_mx=Clamped_Phi(TV(location.x-DX.x,location.y,location.z),phi_nodes),phi_px=Clamped_Phi(TV(location.x+DX.x,location.y,location.z),phi_nodes),
            phi_my=Clamped_Phi(TV(location.x,location.y-DX.y,location.z),phi_nodes),phi_py=Clamped_Phi(TV(location.x,location.y+DX.y,location.z),phi_nodes),
            phi_mz=Clamped_Phi(TV(location.x,location.y,location.z-DX.z),phi_nodes),phi_pz=Clamped_Phi(TV(location.x,location.y,location.z+DX.z),phi_nodes),
            phi_px_py=Clamped_Phi(TV(location.x+DX.x,location.y+DX.y,location.z),phi_nodes),phi_px_my=Clamped_Phi(TV(location.x+DX.x,location.y-DX.y,location.z),phi_nodes),
            phi_mx_py=Clamped_Phi(TV(location.x-DX.x,location.y+DX.y,location.z),phi_nodes),phi_mx_my=Clamped_Phi(TV(location.x-DX.x,location.y-DX.y,location.z),phi_nodes),
            phi_px_pz=Clamped_Phi(TV(location.x+DX.x,location.y,location.z+DX.z),phi_nodes),phi_px_mz=Clamped_Phi(TV(location.x+DX.x,location.y,location.z-DX.z),phi_nodes),
            phi_mx_pz=Clamped_Phi(TV(location.x-DX.x,location.y,location.z+DX.z),phi_nodes),phi_mx_mz=Clamped_Phi(TV(location.x-DX.x,location.y,location.z-DX.z),phi_nodes),
            phi_py_pz=Clamped_Phi(TV(location.x,location.y+DX.y,location.z+DX.z),phi_nodes),phi_py_mz=Clamped_Phi(TV(location.x,location.y+DX.y,location.z-DX.z),phi_nodes),
            phi_my_pz=Clamped_Phi(TV(location.x,location.y-DX.y,location.z+DX.z),phi_nodes),phi_my_mz=Clamped_Phi(TV(location.x,location.y-DX.y,location.z-DX.z),phi_nodes);
        T phix=(phi_px-phi_mx)*one_over_two_dx,phiy=(phi_py-phi_my)*one_over_two_dy,phiz=(phi_pz-phi_mz)*one_over_two_dz,
            phixx=(phi_px-2*phi_curr+phi_mx)*one_over_dx_squared,phiyy=(phi_py-2*phi_curr+phi_my)*one_over_dy_squared,
            phizz=(phi_pz-2*phi_curr+phi_mz)*one_over_dz_squared,phixy=(phi_px_py-phi_px_my-phi_mx_py+phi_mx_my)*one_over_four_dx_dy,
            phixz=(phi_px_pz-phi_px_mz-phi_mx_pz+phi_mx_mz)*one_over_four_dx_dz,phiyz=(phi_py_pz-phi_py_mz-phi_my_pz+phi_my_mz)*one_over_four_dy_dz;
        T denominator=sqrt(sqr(phix)+sqr(phiy)+sqr(phiz));
        T curvature_value;
        if(denominator >= small_number) 
            curvature_value=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx+sqr(phix)*phizz-2*phix*phiz*phixz+sqr(phiz)*phixx+sqr(phiy)*phizz-2*phiy*phiz*phiyz+sqr(phiz)*phiyy)/cube(denominator);
        else curvature_value=LEVELSET_UTILITIES<T>::Sign(phi_curr)*max_curvature;
        return minmag(curvature_value,LEVELSET_UTILITIES<T>::Sign(curvature_value)*max_curvature);}
}
//#####################################################################
template class LEVELSET_OCTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_OCTREE<double>;
#endif
#endif
