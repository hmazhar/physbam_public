//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CHILDREN.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/FACE_LOOKUP_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,2> LEVELSET_QUADTREE<T>::
Hessian(const TV& X,const ARRAY<T>* phi_nodes) const
{   
    TV DX(grid.Leaf_Cell(X)->DX());
    T one_over_dx=1/DX.x,one_over_dy=1/DX.y;T two_phi_center=2*Phi(X,phi_nodes);
    T phi_xx=(Phi(TV(X.x+DX.x,X.y),phi_nodes)-two_phi_center+Phi(TV(X.x-DX.x,X.y),phi_nodes))*sqr(one_over_dx),
        phi_yy=(Phi(TV(X.x,X.y+DX.y),phi_nodes)-two_phi_center+Phi(TV(X.x,X.y-DX.y),phi_nodes))*sqr(one_over_dy),
        phi_xy=(Phi(TV(X.x+DX.x,X.y+DX.y),phi_nodes)-Phi(TV(X.x+DX.x,X.y-DX.y),phi_nodes)-Phi(TV(X.x-DX.x,X.y+DX.y),phi_nodes)+Phi(TV(X.x-DX.x,X.y-DX.y),phi_nodes))
                    *(T).25*one_over_dx*one_over_dy;
    return SYMMETRIC_MATRIX<T,2>(phi_xx,phi_xy,phi_yy);
}   
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,1> LEVELSET_QUADTREE<T>::
Principal_Curvatures(const TV& X,const ARRAY<T>* phi_nodes) const
{
    TV DX(grid.Leaf_Cell(X)->DX());
    TV grad_phi=TV((Phi(TV(X.x+DX.x,X.y),phi_nodes)-Phi(TV(X.x-DX.x,X.y),phi_nodes))/((T)2*DX.x),
                   (Phi(TV(X.x,X.y+DX.y),phi_nodes)-Phi(TV(X.x,X.y-DX.y),phi_nodes))/((T)2*DX.y));
    TV tangent=grad_phi.Perpendicular();T grad_phi_magnitude=tangent.Normalize();
    T curvature=TV::Dot_Product(tangent,Hessian(X)*tangent)/grad_phi_magnitude;
    return VECTOR<T,1>(curvature);
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T> void LEVELSET_QUADTREE<T>::
Compute_Curvature(const T time)
{
    ARRAY<T> phi_ghost(grid.number_of_cells);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);
    ARRAY<T> phi_nodes_ghost(grid.number_of_nodes);LINEAR_INTERPOLATION_DYADIC_HELPER<QUADTREE_GRID<T> >::Interpolate_From_Cells_To_Nodes(grid,phi_ghost,phi_nodes_ghost);
    if(!curvature) curvature=new ARRAY<T>(grid.number_of_nodes);else curvature->Resize(grid.number_of_nodes);

    for(DYADIC_GRID_ITERATOR_CELL<QUADTREE_GRID<T> > iterator(grid,2);iterator.Valid();iterator.Next()){
        TV location=iterator.Location();
        TV DX=iterator.DX();

        T one_over_two_dx=1/(2*DX.x),one_over_two_dy=1/(2*DX.y),one_over_dx_squared=1/sqr(DX.x),one_over_dy_squared=1/sqr(DX.y),one_over_four_dx_dy=1/(4*DX.x*DX.y);
        T max_curvature=1/grid.Minimum_Edge_Length(); // max resolution
        const QUADTREE_CELL<T>* cell=iterator.Cell_Pointer();
        T phi_curr=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y)),
            phi_mx=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y)),
            phi_px=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y)),
            phi_my=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y-DX.y)),
            phi_py=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x,location.y+DX.y)),
            phi_px_py=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y+DX.y)),
            phi_px_my=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x+DX.x,location.y-DX.y)),
            phi_mx_py=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y+DX.y)),
            phi_mx_my=interpolation->From_Close_Cell_Cell(grid,cell,phi_ghost,&phi_nodes_ghost,TV(location.x-DX.x,location.y-DX.y));
        T phix=(phi_px-phi_mx)*one_over_two_dx,phiy=(phi_py-phi_my)*one_over_two_dy,
            phixx=(phi_px-2*phi_curr+phi_mx)*one_over_dx_squared,phiyy=(phi_py-2*phi_curr+phi_my)*one_over_dy_squared,
            phixy=(phi_px_py-phi_px_my-phi_mx_py+phi_mx_my)*one_over_four_dx_dy;
        T denominator=sqrt(sqr(phix)+sqr(phiy));
        T curvature_value;
        if(denominator >= small_number) 
            curvature_value=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx)/cube(denominator);
        else curvature_value=LEVELSET_UTILITIES<T>::Sign(phi_curr)*max_curvature;
        (*curvature)(iterator.Cell_Index())=minmag(curvature_value,LEVELSET_UTILITIES<T>::Sign(curvature_value)*max_curvature);}
}
//#####################################################################
// Function Curvature
//#####################################################################
template<class T> T LEVELSET_QUADTREE<T>::
Curvature(const TV& location,const bool use_precomputed_curvatures_if_possible,const ARRAY<T>* phi_nodes) const
{
    if(curvature && use_precomputed_curvatures_if_possible) return LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,*curvature,location);
    else{ // warning: the following code is incredibly expensive
        TV DX(grid.Leaf_Cell(location)->DX());
        T one_over_two_dx=1/(2*DX.x),one_over_two_dy=1/(2*DX.y),one_over_dx_squared=1/sqr(DX.x),one_over_dy_squared=1/sqr(DX.y),one_over_four_dx_dy=1/(4*DX.x*DX.y);
        T max_curvature=1/grid.Minimum_Edge_Length(); // max resolution

        T phi_curr=Clamped_Phi(TV(location.x,location.y),phi_nodes),
            phi_mx=Clamped_Phi(TV(location.x-DX.x,location.y),phi_nodes),phi_px=Clamped_Phi(TV(location.x+DX.x,location.y),phi_nodes),
            phi_my=Clamped_Phi(TV(location.x,location.y-DX.y),phi_nodes),phi_py=Clamped_Phi(TV(location.x,location.y+DX.y),phi_nodes),
            phi_px_py=Clamped_Phi(TV(location.x+DX.x,location.y+DX.y),phi_nodes),phi_px_my=Clamped_Phi(TV(location.x+DX.x,location.y-DX.y),phi_nodes),
            phi_mx_py=Clamped_Phi(TV(location.x-DX.x,location.y+DX.y),phi_nodes),phi_mx_my=Clamped_Phi(TV(location.x-DX.x,location.y-DX.y),phi_nodes);
        T phix=(phi_px-phi_mx)*one_over_two_dx,phiy=(phi_py-phi_my)*one_over_two_dy,
            phixx=(phi_px-2*phi_curr+phi_mx)*one_over_dx_squared,phiyy=(phi_py-2*phi_curr+phi_my)*one_over_dy_squared,
            phixy=(phi_px_py-phi_px_my-phi_mx_py+phi_mx_my)*one_over_four_dx_dy;
        T denominator=sqrt(sqr(phix)+sqr(phiy));
        T curvature_value;
        if(denominator >= small_number) 
            curvature_value=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx)/cube(denominator);
        else curvature_value=LEVELSET_UTILITIES<T>::Sign(phi_curr)*max_curvature;
        return minmag(curvature_value,LEVELSET_UTILITIES<T>::Sign(curvature_value)*max_curvature);}
}
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class T> struct NEGATIVE_MATERIAL_HELPER{const LEVELSET_QUADTREE<T>* levelset;T* total_area;};
template<class T> static void Negative_Material_Helper(void* data,const QUADTREE_CELL<T>* cell)
{
    NEGATIVE_MATERIAL_HELPER<T>* helper=(NEGATIVE_MATERIAL_HELPER<T>*)data;
    ARRAY<T>& phi=helper->levelset->phi;VECTOR<T,2> dx=cell->DX();
    *helper->total_area+=LEVELSET_UTILITIES<T>::Negative_Cell_Fraction(phi(cell->Node(0)),phi(cell->Node(1)),phi(cell->Node(2)),phi(cell->Node(3)),1)*dx.x*dx.y;
}
template<class T> T LEVELSET_QUADTREE<T>::
Negative_Material() const
{ 
    T total_area=0;
    NEGATIVE_MATERIAL_HELPER<T> helper;helper.levelset=this;helper.total_area=&total_area;
    MAP_QUADTREE_MESH<T>::Map_Cells(grid.uniform_grid,grid.cells,0,&helper,Negative_Material_Helper);
    return total_area;
}
//#####################################################################
// Function Positive_Material
//#####################################################################
 // only works with phi stored at the nodes - but this is always the case for quadtrees
template<class T> T LEVELSET_QUADTREE<T>::
Positive_Material() const
{ 
    RANGE<TV> domain=grid.Domain();
    return domain.Size()-Negative_Material();
}
//#####################################################################
template class LEVELSET_QUADTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_QUADTREE<double>;
#endif
#endif
