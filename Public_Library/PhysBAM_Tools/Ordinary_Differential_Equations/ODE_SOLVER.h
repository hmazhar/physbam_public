//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ODE_SOLVER
//#####################################################################
#ifndef __ODE_SOLVER__
#define __ODE_SOLVER__

#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
namespace PhysBAM{

template<class T,class T2=T>
class ODE_SOLVER
{
public:

    T2 Euler(const NONLINEAR_FUNCTION<T2(T,T2)>& F,const T2 y_not,const T t_not,const T t_final,int n)
    {T2 y=y_not;T h=(t_final-t_not)/n;
    for(int k=1;k<=n;k++){
        T t=t_not+(k-1)*h;
        y+=h*F(t,y);}
    return y;}

    T2 Runge_Kutta_2(const NONLINEAR_FUNCTION<T2(T,T2)>& F,const T2 y_not,const T t_not,const T t_final,int n)
    {T2 y=y_not;T h=(t_final-t_not)/n;
    for(int k=1;k<=n;k++){
        T t=t_not+(k-1)*h;T2 k1=F(t,y),k2=F(t+h,y+h*k1);
        y+=h*(k1+k2)/2;}
    return y;}

    T2 Runge_Kutta_4(const NONLINEAR_FUNCTION<T2(T,T2)>& F,const T2 y_not,const T t_not,const T t_final,int n)
    {T2 y=y_not;T h=(t_final-t_not)/n;
    for(int k=1;k<=n;k++){
        T t=t_not+(k-1)*h;T2 k1=F(t,y),k2=F(t+h/2,y+h*k1/2),k3=F(t+h/2,y+h*k2/2),k4=F(t+h,y+h*k3);
        y+=h*(k1+(T)2*k2+(T)2*k3+k4)/6;}
    return y;}

    T2 Adams_Bashforth(const NONLINEAR_FUNCTION<T2(T,T2)>& F,const T2 y_not,const T t_not,const T t_final,int n)
    {T h=(t_final-t_not)/n;
    T2 ym4=y_not,Fm4=F(t_not,ym4);
    T2 ym3=Runge_Kutta_4(F,ym4,t_not,t_not+h,1),Fm3=F(t_not+h,ym3);
    T2 ym2=Runge_Kutta_4(F,ym3,t_not+h,t_not+2*h,1),Fm2=F(t_not+2*h,ym2);
    T2 ym1=Runge_Kutta_4(F,ym2,t_not+2*h,t_not+3*h,1),Fm1=F(t_not+3*h,ym1);
    T2 y=ym1;
    for(int k=4;k<=n;k++){
        T t=t_not+(k-1)*h;
        y+=h*((T)55*Fm1-(T)59*Fm2+(T)37*Fm3-(T)9*Fm4)/24;
        Fm4=Fm3;Fm3=Fm2;Fm2=Fm1;Fm1=F(t+h,y);}
    return y;}

    T2 Adams_Bashforth_Adams_Moulton(const NONLINEAR_FUNCTION<T2(T,T2)>& F,const T2 y_not,const T t_not,const T t_final,int n)
    {T h=(t_final-t_not)/n;
    T2 ym4=y_not,Fm4=F(t_not,ym4);
    T2 ym3=Runge_Kutta_4(F,ym4,t_not,t_not+h,1),Fm3=F(t_not+h,ym3);
    T2 ym2=Runge_Kutta_4(F,ym3,t_not+h,t_not+2*h,1),Fm2=F(t_not+2*h,ym2);
    T2 ym1=Runge_Kutta_4(F,ym2,t_not+2*h,t_not+3*h,1),Fm1=F(t_not+3*h,ym1);
    T2 y=ym1;
    for(int k=4;k<=n;k++){
        T t=t_not+(k-1)*h;
        T2 y_guess=y+h*((T)55*Fm1-(T)59*Fm2+(T)37*Fm3-(T)9*Fm4)/24,F_guess=F(t+h,y_guess);
        y+=h*((T)9*F_guess+(T)19*Fm1-(T)5*Fm2+Fm3)/24;
        Fm4=Fm3;Fm3=Fm2;Fm2=Fm1;Fm1=F(t+h,y);}
    return y;}

//#####################################################################
};
}
#endif
