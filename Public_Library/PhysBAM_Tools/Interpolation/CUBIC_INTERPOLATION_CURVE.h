//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_INTERPOLATION_CURVE
//#####################################################################
#ifndef __CUBIC_INTERPOLATION_CURVE__
#define __CUBIC_INTERPOLATION_CURVE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,class T2>
class CUBIC_INTERPOLATION_CURVE
{
    struct CONTROL_POINT
    {
        CONTROL_POINT()
        {}

        CONTROL_POINT(T time,T2 value)
        {t=time;x=value;}

        T t;T2 x;
        mutable T2 b,c,d;
    };
    ARRAY<CONTROL_POINT> control_points;
    mutable bool need_compute_coefficients;
public:

    CUBIC_INTERPOLATION_CURVE()
        :need_compute_coefficients(false)
    {}

    int Locate_Interval(const T t) const
    {int lo=0,hi=control_points.m+1;
    while(hi-lo>1){int m=(lo+hi)/2;if(t<control_points(m).t) hi=m;else lo=m;}
    return lo;}

    T2 Value(T t) const
    {if(control_points.m<=1) return control_points.m?control_points(1).x:T2();
    if(need_compute_coefficients) Compute_Coefficients();
    int left_index=Locate_Interval(t);const CONTROL_POINT& pt=control_points(left_index?left_index:1);T h=t-pt.t;
    if(left_index==0 || left_index==control_points.m) return pt.b*h+pt.x;
    return ((pt.d*h+pt.c)*h+pt.b)*h+pt.x;}

    T2 Derivative(T t) const
    {if(control_points.m<=1) return T2();
    if(need_compute_coefficients) Compute_Coefficients();
    int left_index=Locate_Interval(t);const CONTROL_POINT& pt=control_points(left_index?left_index:1);
    if(left_index==0 || left_index==control_points.m) return pt.b;
    T h=t-pt.t;return (pt.d*(3*h)+pt.c*2)*h+pt.b;}

    T2 Second_Derivative(T t) const
    {if(control_points.m<=1) return T2();
    int left_index=Locate_Interval(t);if(left_index==0 || left_index==control_points.m) return T2();
    if(need_compute_coefficients) Compute_Coefficients();
    const CONTROL_POINT& pt=control_points(left_index);return pt.d*((t-pt.t)*6)+pt.c*2;}

    void Add_Control_Point(T t,const T2& value)
    {int index=Locate_Interval(t);need_compute_coefficients=true;
    control_points.Insert(CONTROL_POINT(t,value),index+1);}

    void Compute_Coefficients() const // ACM Algorithm 472
    {need_compute_coefficients=false;
    if(control_points.m<=1) return;T2 r,s=T2();T u=T(),v=T();ARRAY<T> h(control_points.m),h_inv(control_points.m),g(control_points.m);
    for(int i=1;i<control_points.m;i++){const CONTROL_POINT &pt=control_points(i),&ptp1=control_points(i+1);
        h(i)=ptp1.t-pt.t;h_inv(i)=1/h(i);r=(ptp1.x-pt.x)*h_inv(i);pt.c=r-s;s=r;}
    r=s=control_points(1).c=control_points(control_points.m).c=T2();
    for(int i=2;i<control_points.m;i++){const CONTROL_POINT &ptm1=control_points(i-1),&pt=control_points(i),&ptp1=control_points(i+1);
        pt.c+=u*ptm1.c;g(i)=1/((ptm1.t-ptp1.t)*2-u*v);v=h(i);u=v*g(i);}
    for(int i=control_points.m-1;i>1;i--){const CONTROL_POINT &pt=control_points(i),&ptp1=control_points(i+1);
        pt.c=(h(i)*ptp1.c-pt.c)*g(i);}
    for(int i=1;i<control_points.m;i++){const CONTROL_POINT &pt=control_points(i),&ptp1=control_points(i+1);
        pt.b=(ptp1.x-pt.x)*h_inv(i)-(pt.c*2+ptp1.c)*h(i);pt.d=(ptp1.c-pt.c)*h_inv(i);pt.c*=3;}
    const CONTROL_POINT &pt=control_points(control_points.m),&ptm1=control_points(control_points.m-1);
    pt.b=(ptm1.d*(3*h(control_points.m-1))+ptm1.c*2)*h(control_points.m-1)+ptm1.b;
}
//#####################################################################
};
}
#endif
