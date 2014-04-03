//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_HAMILTON_JACOBI_ENO  
//##################################################################### 
#ifndef __ADVECTION_HAMILTON_JACOBI_ENO__
#define __ADVECTION_HAMILTON_JACOBI_ENO__

#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class ADVECTION_HAMILTON_JACOBI_ENO:public ADVECTION_SEPARABLE_UNIFORM<T_GRID,T2>
{
private:
    typedef typename T_GRID::SCALAR T;

    int order;
public:

    ADVECTION_HAMILTON_JACOBI_ENO()
    {
        Set_Order();
    }

    void Set_Order(const int order_input=3)  PHYSBAM_OVERRIDE
    {order=order_input;assert(order >=1 && order <=3);}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) PHYSBAM_OVERRIDE;
    void Advection_Solver(const int m_start,const int m_end,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx);
//#####################################################################
};
}    
#endif
