//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO  
//##################################################################### 
#ifndef __ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO__
#define __ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO__

#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class T2>
class ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO:public ADVECTION_SEPARABLE_UNIFORM<T_input,T2>
{
    typedef T_input T;

    bool compute_epsilon;
    T epsilon;
public:

    ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO()
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
template<class T,class T2> void ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO<T,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    T one_over_dx=1/dx;
    ARRAY<T2,VECTOR<int,1> > D1(-2,m+3); // 1st divided difference
    for(int i=-2;i<=m+3;i++) D1(i)=Z(i);

    if(compute_epsilon){
        epsilon=1e-6*sqr(maxabs(D1)); 
        if(epsilon == 0) epsilon=1e-6;} // epsilon=0 implies all v_i=0 and all u_Zx=0

    for(int i=1;i<=m;i++){
        if(u(i) > 0){
            T2 flux_left=WENO(D1(i-3),D1(i-2),D1(i-1),D1(i),D1(i+1),epsilon),flux_right=WENO(D1(i-2),D1(i-1),D1(i),D1(i+1),D1(i+2),epsilon);
            u_Zx(i)=u(i)*(flux_right-flux_left)*one_over_dx;}
        else{ 
            T2 flux_left=WENO(D1(i+2),D1(i+1),D1(i),D1(i-1),D1(i-2),epsilon),flux_right=WENO(D1(i+3),D1(i+2),D1(i+1),D1(i),D1(i-1),epsilon);
            u_Zx(i)=u(i)*(flux_right-flux_left)*one_over_dx;}}
}
//#####################################################################
}
#endif

