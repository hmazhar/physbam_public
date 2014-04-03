//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTEGRATION
//#####################################################################
#ifndef __INTEGRATION__
#define __INTEGRATION__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <cmath>
namespace PhysBAM{

using ::std::sqrt;

template<class T>
class INTEGRATION
{
public:

    T Midpoint(const NONLINEAR_FUNCTION<T(T)>& F,const T a,const T b,int n)
    {T h=(b-a)/n,sum=0;
    for(int i=1;i<=n;i++) sum+=F(a+(i-.5)*h);
    return h*sum;}

    T Trapezoid(const NONLINEAR_FUNCTION<T(T)>& F,const T a,const T b,int n)
    {T h=(b-a)/n,sum=0,F_left=F(a);
    for(int i=1;i<=n;i++){T F_right=F(a+i*h);sum+=(F_left+F_right);F_left=F_right;}
    return h/2*sum;}

    T Simpson(const NONLINEAR_FUNCTION<T(T)>& F,const T a,const T b,int n)
    {T h=(b-a)/n,sum=0,F_left=F(a);
    for(int i=1;i<=n;i++){T x=a+(i-1)*h,F_right=F(x+h);sum+=(F_left+4*F(x+h/2)+F_right);F_left=F_right;}
    return h/6*sum;}

    T Gaussian(const NONLINEAR_FUNCTION<T(T)>& F,const T a,const T b,int n)
    {T h=(b-a)/n,c=1/(2*sqrt(3)),sum=0;
    for(int i=1;i<=n;i++){T m=a+(i-.5)*h,d=c*h;sum+=(F(m-d)+F(m+d));}
    return h/2*sum;}

//#####################################################################
};
}
#endif
