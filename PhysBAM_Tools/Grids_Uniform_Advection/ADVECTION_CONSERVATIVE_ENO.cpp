//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_CONSERVATIVE_ENO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_input,class T2> ADVECTION_CONSERVATIVE_ENO<T_input,T2>::
ADVECTION_CONSERVATIVE_ENO()
{
    Set_Order();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_input,class T2> ADVECTION_CONSERVATIVE_ENO<T_input,T2>::
~ADVECTION_CONSERVATIVE_ENO()
{
}
//#####################################################################
// Function Set_Order
//#####################################################################
template<class T_input,class T2> void ADVECTION_CONSERVATIVE_ENO<T_input,T2>::
Set_Order(const int order_input)
{
    order=order_input;
    assert(order>=1 && order<=3);
}
//#####################################################################
// Function Advection_Solver
//#####################################################################
// finds (uZ)_x with Local Lax Friedrichs
template<class T,class T2> void ADVECTION_CONSERVATIVE_ENO<T,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    T one_over_dx=1/dx,one_over_two_dx=.5*one_over_dx,one_over_three_dx=one_third*one_over_dx;
    ARRAY<T2,VECTOR<int,1> > DZ1(-2,m+3),DZ2(-2,m+2),DZ3(-2,m+1); // divided differences
    for(int i=-2;i<=m+3;i++) DZ1(i)=Z(i); 
    if(order >= 2) for(int i=-2;i<=m+2;i++) DZ2(i)=(DZ1(i+1)-DZ1(i))*one_over_two_dx;     
    if(order == 3) for(int i=-2;i<=m+1;i++) DZ3(i)=(DZ2(i+1)-DZ2(i))*one_over_three_dx;
    ARRAY<T2,VECTOR<int,1> > DUZ1(-2,m+3),DUZ2(-2,m+2),DUZ3(-2,m+1); // divided differences
    for(int i=-2;i<=m+3;i++) DUZ1(i)=u(i)*Z(i); 
    if(order >= 2) for(int i=-2;i<=m+2;i++) DUZ2(i)=(DUZ1(i+1)-DUZ1(i))*one_over_two_dx;     
    if(order == 3) for(int i=-2;i<=m+1;i++) DUZ3(i)=(DUZ2(i+1)-DUZ2(i))*one_over_three_dx;

    ARRAY<T2,VECTOR<int,1> > flux(0,m); // flux is to the right of each point 
    if(order == 1) for(int i=0;i<=m;i++){
        T2 alpha=maxabs(u(i),u(i+1));
        T2 flux_left=DUZ1(i)+alpha*DZ1(i);
        T2 flux_right=DUZ1(i+1)-alpha*DZ1(i+1);
        flux(i)=.5*(flux_left+flux_right);}
    else if(order == 2) for(int i=0;i<=m;i++){
        T2 alpha=maxabs(u(i),u(i+1));
        T2 flux_left=ENO(dx,DUZ1(i)+alpha*DZ1(i),DUZ2(i-1)+alpha*DZ2(i-1),DUZ2(i)+alpha*DZ2(i));
        T2 flux_right=ENO(dx,DUZ1(i+1)-alpha*DZ1(i+1),-(DUZ2(i+1)-alpha*DZ2(i+1)),-(DUZ2(i)-alpha*DZ2(i)));
        flux(i)=.5*(flux_left+flux_right);}
    else if(order == 3) for(int i=0;i<=m;i++){
        T2 alpha=maxabs(u(i),u(i+1));
        T2 flux_left=ENO(dx,DUZ1(i)+alpha*DZ1(i),DUZ2(i-1)+alpha*DZ2(i-1),DUZ2(i)+alpha*DZ2(i),DUZ3(i-2)+alpha*DZ3(i-2),DUZ3(i-1)+alpha*DZ3(i-1),DUZ3(i)+alpha*DZ3(i));
        T2 flux_right=ENO(dx,DUZ1(i+1)-alpha*DZ1(i+1),-(DUZ2(i+1)-alpha*DZ2(i+1)),-(DUZ2(i)-alpha*DZ2(i)),DUZ3(i+1)-alpha*DZ3(i+1),DUZ3(i)-alpha*DZ3(i),DUZ3(i-1)-alpha*DZ3(i-1));
        flux(i)=.5*(flux_left+flux_right);}

    for(int i=1;i<=m;i++) u_Zx(i)=(flux(i)-flux(i-1))*one_over_dx;
}
//#####################################################################
template class ADVECTION_CONSERVATIVE_ENO<GRID<VECTOR<float,1> >,float>;
template class ADVECTION_CONSERVATIVE_ENO<GRID<VECTOR<float,2> >,float>;
template class ADVECTION_CONSERVATIVE_ENO<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_CONSERVATIVE_ENO<GRID<VECTOR<double,1> >,double>;
template class ADVECTION_CONSERVATIVE_ENO<GRID<VECTOR<double,2> >,double>;
template class ADVECTION_CONSERVATIVE_ENO<GRID<VECTOR<double,3> >,double>;
#endif
