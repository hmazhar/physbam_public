//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYNOMIAL
//#####################################################################
#ifndef __POLYNOMIAL__
#define __POLYNOMIAL__

#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <cassert>
namespace PhysBAM{

template<class T>
class POLYNOMIAL:public NONLINEAR_FUNCTION<T(T)>
{
public:
    int degree;
    VECTOR_ND<T> c; // coefficients

    POLYNOMIAL(const int degree_input)
        :degree(degree_input),c(degree+1)
    {}

    T operator()(const T x_input) const PHYSBAM_OVERRIDE
    {T y=c(degree+1);for(int i=degree;i>=1;i--) y=c(i)+x_input*y;return y;}

    void Compute_Coefficients(const VECTOR_ND<T>& x,const VECTOR_ND<T>& y)
    {assert(x.n == degree+1 && y.n == degree+1);MATRIX_MXN<T> A(degree+1,degree+1);
    for(int i=1;i<=degree+1;i++){A(i,1)=1;for(int j=2;j<=degree+1;j++)A(i,j)=x(i)*A(i,j-1);}
    VECTOR_ND<T> c_temp=A.PLU_Solve(y);for(int i=1;i<=degree+1;i++)c(i)=c_temp(i);}

//#####################################################################
};
}
#endif
