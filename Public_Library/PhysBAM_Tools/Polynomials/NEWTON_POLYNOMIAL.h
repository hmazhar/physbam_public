//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEWTON_POLYNOMIAL
//#####################################################################
#ifndef __NEWTON_POLYNOMIAL__
#define __NEWTON_POLYNOMIAL__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T>
class NEWTON_POLYNOMIAL:public NONLINEAR_FUNCTION<T(T)>
{
public:
    int degree;
    VECTOR_ND<T> c; // coefficients
    VECTOR_ND<T> x; // data points

    NEWTON_POLYNOMIAL(const int degree_input)
        :degree(degree_input),c(degree+1),x(degree+1)
    {}

    T operator()(const T x_input) const PHYSBAM_OVERRIDE
    {T y=c(degree+1);for(int i=degree;i>=1;i--) y=c(i)+(x_input-x(i))*y;return y;}

    void Compute_Coefficients(const VECTOR_ND<T>& x_input,const VECTOR_ND<T>& y_input)
    {int i,j;assert(x_input.n == degree+1 && y_input.n == degree+1);for(i=1;i<=degree+1;i++) x(i)=x_input(i);
    VECTOR_ND<T> f(degree+1);for(i=1;i<=degree+1;i++) f(i)=y_input(i); c(1)=f(1);
    for(j=1;j<=degree;j++){for(i=1;i<=degree+1-j;i++) f(i)=(f(i+1)-f(i))/(x(i+1+j-1)-x(i)); c(j+1)=f(1);}}

//#####################################################################
};
}
#endif
