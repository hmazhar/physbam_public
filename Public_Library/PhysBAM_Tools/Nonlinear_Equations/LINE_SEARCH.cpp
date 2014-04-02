//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <cmath>
using std::abs;
using namespace PhysBAM;
//#####################################################################
// Function Update_Interval
//#####################################################################
template<class T> void LINE_SEARCH<T>::
Update_Interval(BRACKET& s,T d,T Fd) const
{
    if(s.m>d){exchange(s.m,d);exchange(s.Fm,Fd);} // a < m < d < b
    T am=min(s.Fa,s.Fm),bd=min(s.Fb,Fd);
    if(am<bd){s.b=d;s.Fb=Fd;} // a m d
    else{s.a=s.m;s.m=d;s.Fa=s.Fm;s.Fm=Fd;} // m d b
}
//#####################################################################
// Function Compute_Quadratic_Minimum
//#####################################################################
template<class T> int LINE_SEARCH<T>::
Compute_Quadratic_Minimum(const BRACKET& s,T& x,T tolerance) const // 0=terminate, 1=fail, 2=good
{
    T m=(s.Fa-s.Fm)*(s.m-s.b),n=(s.Fm-s.Fb)*(s.a-s.m),p=2*(n-m);
    if(p<0) return 1;
    if(p<=tolerance*(abs(s.Fa)+abs(s.Fm)+abs(s.Fb))*(abs(s.a)+abs(s.b))) return 0;
    T r=(s.a+s.m)*n-(s.m+s.b)*m;
    x=r/p;
    if(x<=s.a || x>=s.m || x==s.m) return 1;
    return 2;
}
//#####################################################################
// Function Best_Value
//#####################################################################
template<class T> T LINE_SEARCH<T>::
Best_Value(const BRACKET& s) const
{
    if(s.Fa<=s.Fm) return (s.Fa<=s.Fb)?s.a:s.b;
    return (s.Fm<=s.Fb)?s.m:s.b;
}
//#####################################################################
// Function Line_Search_Quadratic_Golden_Section
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Quadratic_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance,T quadratic_tolerance) const
{
    T m=(T).5*(a+b),t;
    BRACKET s={a,m,b,F(a),F(m),F(b)};
    int i;
    for(i=1;i<=max_iterations;i++){
        if(abs(s.b-s.a)<interval_tolerance*(abs(s.a)+abs(s.b))) break;
        int r=Compute_Quadratic_Minimum(s,t,quadratic_tolerance);
        if(!r) break;
        if(r==1){
            if(s.m-s.a>s.b-s.m) t=(T).5*(s.a+s.m);
            else t=(T).5*(s.m+s.b);}
        Update_Interval(s,t,F(t));}
    x=Best_Value(s);
    return i<=max_iterations;
}
//#####################################################################
// Function Line_Search_Quadratic_Golden_Section
//#####################################################################
template<class T> bool LINE_SEARCH<T>::
Line_Search_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance) const
{
    PHYSBAM_ASSERT(a<=b);
    T tau=(T).5*(sqrt((T)5)-1),m=a+tau*(b-a),t;
    BRACKET s={a,m,b,F(a),F(m),F(b)};
    int i;
    interval_tolerance*=abs(s.a)+abs(s.b);
    for(i=1;i<=max_iterations;i++){
        if(abs(s.b-s.a)<interval_tolerance) break;
        t=s.a+s.b-s.m;
        Update_Interval(s,t,F(t));}
    x=Best_Value(s);
    return i<=max_iterations;
}
template class LINE_SEARCH<float>;
template class LINE_SEARCH<double>;
