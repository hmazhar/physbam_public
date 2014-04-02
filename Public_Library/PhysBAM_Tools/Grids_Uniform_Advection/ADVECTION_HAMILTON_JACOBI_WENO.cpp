//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_HAMILTON_JACOBI_WENO  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_WENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
using namespace PhysBAM;
//##################################################################### 
// Function Advection_Solver
//#####################################################################
template<class T_GRID,class T2> void ADVECTION_HAMILTON_JACOBI_WENO<T_GRID,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{      
    int i;T one_over_dx=1/dx;
    ARRAY<T2,VECTOR<int,1> > D1(-2,m+2);for(i=-2;i<=m+2;i++) D1(i)=(Z(i+1)-Z(i))*one_over_dx;

    if(compute_epsilon){
        epsilon=(T)1e-6*sqr(D1.Maxabs());
        if(epsilon == 0) epsilon=(T)1e-6;} // epsilon=0 implies all v_i=0 and all u_Zx=0  

    for(i=1;i<=m;i++){
        if(u(i) > 0) u_Zx(i)=u(i)*WENO(D1(i-3),D1(i-2),D1(i-1),D1(i),D1(i+1),epsilon);
        else u_Zx(i)=u(i)*WENO(D1(i+2),D1(i+1),D1(i),D1(i-1),D1(i-2),epsilon);}
}
//#####################################################################
template class ADVECTION_HAMILTON_JACOBI_WENO<GRID<VECTOR<float,1> >,float>;
template class ADVECTION_HAMILTON_JACOBI_WENO<GRID<VECTOR<float,2> >,float>;
template class ADVECTION_HAMILTON_JACOBI_WENO<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_HAMILTON_JACOBI_WENO<GRID<VECTOR<double,1> >,double>;
template class ADVECTION_HAMILTON_JACOBI_WENO<GRID<VECTOR<double,2> >,double>;
template class ADVECTION_HAMILTON_JACOBI_WENO<GRID<VECTOR<double,3> >,double>;
#endif
