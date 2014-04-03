//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_CENTRAL
//##################################################################### 
#ifndef __ADVECTION_CENTRAL__
#define __ADVECTION_CENTRAL__

#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class ADVECTION_CENTRAL:public ADVECTION_SEPARABLE_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;
public:
    ADVECTION_CENTRAL()
    {}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Advection_Solver
//#####################################################################
template<class T_GRID,class T2> void ADVECTION_CENTRAL<T_GRID,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    T one_over_two_dx=1/(2*dx);
    for(int i=1;i<=m;i++) u_Zx(i)=u(i)*(Z(i+1)-Z(i-1))*one_over_two_dx;
}
//#####################################################################
}    
#endif
