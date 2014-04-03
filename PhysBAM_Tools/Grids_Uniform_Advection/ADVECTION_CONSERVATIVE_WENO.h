//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_CONSERVATIVE_WENO  
//##################################################################### 
#ifndef __ADVECTION_CONSERVATIVE_WENO__
#define __ADVECTION_CONSERVATIVE_WENO__

#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class T2>
class ADVECTION_CONSERVATIVE_WENO:public ADVECTION_SEPARABLE_UNIFORM<T_input,T2>
{
    typedef T_input T;

    bool compute_epsilon;
    double epsilon;
public:

    ADVECTION_CONSERVATIVE_WENO()
    {
        Compute_Epsilon();
        Set_Epsilon();
    }

    void Set_Epsilon(const T epsilon_input=1e-6) PHYSBAM_OVERRIDE
    {compute_epsilon=false;epsilon=epsilon_input;}
    
    void Compute_Epsilon() PHYSBAM_OVERRIDE
    {compute_epsilon=true;}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) PHYSBAM_OVERRIDE;
//#####################################################################
};   
//#####################################################################
// Function Advection_Solver
//#####################################################################
// finds (uZ)_x with Local Lax Friedrichs
template<class T,class T2> void ADVECTION_CONSERVATIVE_WENO<T,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    ARRAY<T2,VECTOR<int,1> > DZ1(-2,m+3); // 1st divided difference
    for(int i=-2;i<=m+3;i++) DZ1(i)=Z(i);
    ARRAY<T2,VECTOR<int,1> > DUZ1(-2,m+3); // 1st divided difference
    for(int i=-2;i<=m+3;i++) DUZ1(i)=u(i)*Z(i);

    if(compute_epsilon){
        epsilon=1e-6*sqr(maxabs(DUZ1)); // only DUZ used to find epsilon 
        if(epsilon == 0) epsilon=1e-6;} // epsilon=0 implies all v_i=0 and all u_Zx=0

    ARRAY<T2,VECTOR<int,1> > flux(0,m); // flux is to the right of each point 
    for(int i=0;i<=m;i++){
        T2 alpha=maxabs(u(i),u(i+1));           
        T2 fluxleft=WENO(DUZ1(i-2)+alpha*DZ1(i-2),DUZ1(i-1)+alpha*DZ1(i-1),DUZ1(i)+alpha*DZ1(i),DUZ1(i+1)+alpha*DZ1(i+1),DUZ1(i+2)+alpha*DZ1(i+2),epsilon);
        T2 fluxright=WENO(DUZ1(i+3)-alpha*DZ1(i+3),DUZ1(i+2)-alpha*DZ1(i+2),DUZ1(i+1)-alpha*DZ1(i+1),DUZ1(i)-alpha*DZ1(i),DUZ1(i-1)-alpha*DZ1(i-1),epsilon);
        flux(i)=.5*(fluxleft+fluxright);}

    T one_over_dx=1/dx;
    for(int i=1;i<=m;i++) u_Zx(i)=(flux(i)-flux(i-1))*one_over_dx;
}
//#####################################################################
}
#endif

