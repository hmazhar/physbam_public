//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGRANGE_POLYNOMIAL
//#####################################################################
#ifndef __LAGRANGE_POLYNOMIAL__
#define __LAGRANGE_POLYNOMIAL__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T>
class LAGRANGE_POLYNOMIAL:public NONLINEAR_FUNCTION<T(T)>
{
public:
    int degree;
    VECTOR_ND<T> c; // coefficients
    VECTOR_ND<T> x; // data points
    VECTOR_ND<T> d; // denominators

    LAGRANGE_POLYNOMIAL(const int degree_input)
        :degree(degree_input),c(degree+1),x(degree+1),d(degree+1)
    {}

    T operator()(const T x_input) const PHYSBAM_OVERRIDE
    {T result=0;
    for(int i=1;i<=degree+1;i++){
        T phi=1;for(int j=1;j<=degree+1;j++) if(i != j) phi*=(x_input-x(j));
        result+=c(i)*phi/d(i);}
    return result;}

    void Compute_Coefficients(const VECTOR_ND<T>& x_input,const VECTOR_ND<T>& y_input)
    {int i,j;assert(x_input.n == degree+1 && y_input.n == degree+1);
    for(i=1;i<=degree+1;i++){x(i)=x_input(i);c(i)=y_input(i);}
    for(i=1;i<=degree+1;i++){d(i)=1;for(j=1;j<=degree+1;j++) if(i != j) d(i)*=(x(i)-x(j));}}

//#####################################################################
};
}
#endif
