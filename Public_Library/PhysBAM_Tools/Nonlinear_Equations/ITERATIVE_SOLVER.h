//#####################################################################
// Copyright 2003-2008, Ronald Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ITERATIVE_SOLVER
//##################################################################### 
#ifndef __ITERATIVE_SOLVER__
#define __ITERATIVE_SOLVER__

namespace PhysBAM{

template<class F> class NONLINEAR_FUNCTION;
template<class T> struct PARAMETER_SPACE;

template<class T>
class ITERATIVE_SOLVER 
{
public:
    int iterations,max_iterations;
    T tolerance;

    ITERATIVE_SOLVER()
        :iterations(0),max_iterations(1000000),tolerance((T)1e-14)
    {}
    
//#####################################################################
    T Bisection_Root(NONLINEAR_FUNCTION<T(T)>& F,T a,T b);
    T Newton_Root(NONLINEAR_FUNCTION<T(T)>& F,T x0);
    T Secant_Root(NONLINEAR_FUNCTION<T(T)>& F,T x0,T x1);
    T Bisection_Secant_Root(NONLINEAR_FUNCTION<T(T)>& F,T a,T b);
    T Bisection_Secant_Root_For_Thin_Shells(NONLINEAR_FUNCTION<T(T)>& F,T a,T b);
    T Bisection_Newton_Root(NONLINEAR_FUNCTION<T(T)>& F,T a,T b);
    T Golden_Minimum(NONLINEAR_FUNCTION<T(T)>& F,T a,T b);
    T Parabolic_Minimum(NONLINEAR_FUNCTION<T(T)>& F,T x0,T x1,T x2);
    T Golden_Parabolic_Minimum(NONLINEAR_FUNCTION<T(T)>& F,T a,T b);
    void Steepest_Decent(NONLINEAR_FUNCTION<T(T,T)>& F,T& x,T& y,const T alpha_max);
    void Conjugate_Gradient(NONLINEAR_FUNCTION<T(T,T)>& F,T& x,T& y,const T alpha_max);
    void Conjugate_Gradient(NONLINEAR_FUNCTION<T(PARAMETER_SPACE<T>)>& F,PARAMETER_SPACE<T>& x,const T alpha_max,const int restart_iterations);
//#####################################################################
};   
}
#endif
