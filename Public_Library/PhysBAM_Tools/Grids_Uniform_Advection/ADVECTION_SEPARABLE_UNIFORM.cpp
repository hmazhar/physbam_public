//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEPARABLE_UNIFORM  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
using namespace PhysBAM;
//#####################################################################
// Function Update_Advection_Equation_Helper
//#####################################################################
template<class T,class T2,class T_AVERAGING> void
Update_Advection_Equation_Helper(ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T2,T_AVERAGING>& advection,const GRID<VECTOR<T,1> >& grid,ARRAYS_ND_BASE<VECTOR<T2,1> >& Z,
    const ARRAYS_ND_BASE<VECTOR<T2,1> >& Z_ghost,const ARRAYS_ND_BASE<VECTOR<VECTOR<T,1>,1> >& V,const T dt,const T time)
{       
    int i;int m=grid.counts.x;T dx=grid.dX.x;ARRAY<T2,VECTOR<int,1> > rhs(1,m);

    ARRAY<T2,VECTOR<int,1> > Z_1d_x(-2,m+3);ARRAY<T,VECTOR<int,1> > u_1d(V.domain.min_corner.x,V.domain.max_corner.x);
    for(i=-2;i<=m+3;i++) Z_1d_x(i)=Z_ghost(i);
    for(i=V.domain.min_corner.x;i<=V.domain.max_corner.x;i++) u_1d(i)=V(i).x;
    advection.Advection_Solver(m,dx,Z_1d_x,u_1d,rhs);

    for(i=1;i<=m;i++) Z(i)-=dt*rhs(i);
}
//#####################################################################
// Function Update_Advection_Equation_Helper
//#####################################################################
template<class T,class T2,class T_AVERAGING> void
Update_Advection_Equation_Helper(ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,2> >,T2,T_AVERAGING>& advection,const GRID<VECTOR<T,2> >& grid,ARRAYS_ND_BASE<VECTOR<T2,2> >& Z,
    const ARRAYS_ND_BASE<VECTOR<T2,2> >& Z_ghost,const ARRAYS_ND_BASE<VECTOR<VECTOR<T,2>,2> >& V,const T dt,const T time)
{
    int i,j;int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;ARRAY<T2,VECTOR<int,2> > rhs(1,m,1,n);

    {ARRAY<T2,VECTOR<int,1> > Z_1d_x(-2,m+3),u_Zx_1d(1,m);ARRAY<T,VECTOR<int,1> > u_1d(V.domain.min_corner.x,V.domain.max_corner.x);
    for(j=1;j<=n;j++){
        for(i=-2;i<=m+3;i++) Z_1d_x(i)=Z_ghost(i,j);
        for(i=V.domain.min_corner.x;i<=V.domain.max_corner.x;i++) u_1d(i)=V(i,j).x;
        advection.Advection_Solver(m,dx,Z_1d_x,u_1d,u_Zx_1d);
        for(i=1;i<=m;i++) rhs(i,j)=u_Zx_1d(i);}}

    {ARRAY<T2,VECTOR<int,1> > Z_1d_y(-2,n+3),v_Zy_1d(1,n);ARRAY<T,VECTOR<int,1> > v_1d(V.domain.min_corner.y,V.domain.max_corner.y);
    for(i=1;i<=m;i++){
        for(j=-2;j<=n+3;j++) Z_1d_y(j)=Z_ghost(i,j);
        for(j=V.domain.min_corner.y;j<=V.domain.max_corner.y;j++) v_1d(j)=V(i,j).y;
        advection.Advection_Solver(n,dy,Z_1d_y,v_1d,v_Zy_1d);
        for(j=1;j<=n;j++) rhs(i,j)+=v_Zy_1d(j);}}

    for(i=1;i<=m;i++) for(j=1;j<=n;j++) Z(i,j)-=dt*rhs(i,j); 
}
//#####################################################################
// Function Update_Advection_Equation_Helper
//#####################################################################
template<class T,class T2,class T_AVERAGING> void
Update_Advection_Equation_Helper(ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,3> >,T2,T_AVERAGING>& advection,const GRID<VECTOR<T,3> >& grid,ARRAYS_ND_BASE<VECTOR<T2,3> >& Z,
    const ARRAYS_ND_BASE<VECTOR<T2,3> >& Z_ghost,const ARRAYS_ND_BASE<VECTOR<VECTOR<T,3>,3> >& V,const T dt,const T time)
{       
    int i,j,ij;int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;ARRAY<T2,VECTOR<int,3> > rhs(1,m,1,n,1,mn);
    
    {ARRAY<T2,VECTOR<int,1> > Z_1d_x(-2,m+3),u_Zx_1d(1,m);ARRAY<T,VECTOR<int,1> > u_1d(V.domain.min_corner.x,V.domain.max_corner.x);
    for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++){
        for(i=-2;i<=m+3;i++) Z_1d_x(i)=Z_ghost(i,j,ij);
        for(i=V.domain.min_corner.x;i<=V.domain.max_corner.x;i++) u_1d(i)=V(i,j,ij).x;
        advection.Advection_Solver(m,dx,Z_1d_x,u_1d,u_Zx_1d);
        for(i=1;i<=m;i++) rhs(i,j,ij)=u_Zx_1d(i);}}

    {ARRAY<T2,VECTOR<int,1> > Z_1d_y(-2,n+3),v_Zy_1d(1,n);ARRAY<T,VECTOR<int,1> > v_1d(V.domain.min_corner.y,V.domain.max_corner.y);
    for(i=1;i<=m;i++) for(ij=1;ij<=mn;ij++){
        for(j=-2;j<=n+3;j++) Z_1d_y(j)=Z_ghost(i,j,ij);
        for(j=V.domain.min_corner.y;j<=V.domain.max_corner.y;j++) v_1d(j)=V(i,j,ij).y;
        advection.Advection_Solver(n,dy,Z_1d_y,v_1d,v_Zy_1d);
        for(j=1;j<=n;j++) rhs(i,j,ij)+=v_Zy_1d(j);}}

    {ARRAY<T2,VECTOR<int,1> > Z_1d_z(-2,mn+3),w_Zz_1d(1,mn);ARRAY<T,VECTOR<int,1> > w_1d(V.domain.min_corner.z,V.domain.max_corner.z);
    for(i=1;i<=m;i++) for(j=1;j<=n;j++){
        for(ij=-2;ij<=mn+3;ij++) Z_1d_z(ij)=Z_ghost(i,j,ij);
        for(ij=V.domain.min_corner.z;ij<=V.domain.max_corner.z;ij++) w_1d(ij)=V(i,j,ij).z;
        advection.Advection_Solver(mn,dz,Z_1d_z,w_1d,w_Zz_1d);
        for(ij=1;ij<=mn;ij++) rhs(i,j,ij)+=w_Zz_1d(ij);}}
        
    for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++) Z(i,j,ij)-=dt*rhs(i,j,ij); 
}
//#####################################################################
// Function Update_Advection_Equation_Node
//#####################################################################
template<class T_GRID,class T2,class T_AVERAGING> void ADVECTION_SEPARABLE_UNIFORM<T_GRID,T2,T_AVERAGING>::
Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    assert(!Z_min && !Z_max);
    Update_Advection_Equation_Helper(*this,grid,Z,Z_ghost,V,dt,time);
}
//#####################################################################
// Function Update_Advection_Equation_Cell
//#####################################################################
template<class T_GRID,class T2,class T_AVERAGING> void ADVECTION_SEPARABLE_UNIFORM<T_GRID,T2,T_AVERAGING>::
Update_Advection_Equation_Cell(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    assert(!Z_min && !Z_max);
    Update_Advection_Equation_Helper(*this,grid,Z,Z_ghost,V,dt,time);
}
//#####################################################################
// Function Update_Advection_Equation_Cell
//#####################################################################
template<class T_GRID,class T2,class T_AVERAGING> void ADVECTION_SEPARABLE_UNIFORM<T_GRID,T2,T_AVERAGING>::
Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,const T_FACE_LOOKUP& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    assert(!Z_min && !Z_max);
    
    T_ARRAYS_VECTOR V_cell(grid.Domain_Indices(3));T_AVERAGING averaging;
    
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) V_cell(iterator.Cell_Index())=averaging.Face_To_Cell_Vector(grid,iterator.Cell_Index(),V);
    
    Update_Advection_Equation_Helper(*this,grid,Z,Z_ghost,V_cell,dt,time);
}
//#####################################################################
// Function Update_Advection_Equation_Cell
//#####################################################################
template<class T_GRID,class T2,class T_AVERAGING> void ADVECTION_SEPARABLE_UNIFORM<T_GRID,T2,T_AVERAGING>::
Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,const T_FACE_LOOKUP& V,T_BOUNDARY& boundary,
    const T dt,const T time,const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    assert(!Z_min && !Z_max);
    
    ARRAY<TV,FACE_INDEX<TV::m> > V_face(grid,0);T_AVERAGING averaging;

    for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(grid);iterator.Valid();iterator.Next()) V_face(iterator.Full_Index())=averaging.Face_To_Face_Vector(grid,iterator.Full_Index(),V);

    for(int i=1;i<=TV::m;i++){
        GRID<TV> node_grid(grid.Get_Axis_X_Face_Grid(i));
        Update_Advection_Equation_Helper(*this,node_grid,Z.Component(i),Z_ghost.V_face.Component(i),V_face.Component(i),dt,time);}
}
//#####################################################################
template class ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<float,1> >,float>;
template class ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<float,2> >,float>;
template class ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<double,1> >,double>;
template class ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<double,2> >,double>;
template class ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<double,3> >,double>;
#endif
