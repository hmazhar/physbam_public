//#####################################################################
// Copyright 2002, 2003, Zhaosheng Bao, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/exchange_sort.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <cassert>
#include <cmath>
using namespace PhysBAM;

using ::std::abs;
using ::std::pow;

//#####################################################################
// Constructor
//#####################################################################
template<class T> CUBIC<T>::
CUBIC(const T c3_input,const T c2_input,const T c1_input,const T c0_input)
    :c3(c3_input),c2(c2_input),c1(c1_input),c0(c0_input),roots(0),error_tolerance((T)1e-14),number_of_extrema(0)
{}
//#####################################################################
// Function operator()
//#####################################################################
template<class T> T CUBIC<T>::
operator()(const T x) const PHYSBAM_OVERRIDE
{
    return Value(x);
}
//#####################################################################
// Function Prime
//#####################################################################
template<class T> T CUBIC<T>::
Prime(const T x) const PHYSBAM_OVERRIDE
{
    return (3*c3*x+2*c2)*x+c1;
}
//#####################################################################
// Function Compute_Roots_Noniterative_In_Interval
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots_Noniterative_In_Interval(const T xmin,const T xmax)
{
    Compute_Roots_Noniterative();
    if(roots == 1){if(root1 < xmin || root1 > xmax) roots=0;}
    else{ // roots=3
        if(root3 < xmin || root1 > xmax){roots=0;return;}
        if(root2 < xmin){if(root3 > xmax){roots=0;return;}else{roots=1;root1=root3;return;}}
        if(root2 > xmax){if(root1 < xmin){roots=0;return;}else{roots=1;return;}}
        if(root3 < xmax){if(root1 > xmin) return;else{roots=2;root1=root2;root2=root3;}}
        if(root1 > xmin){roots=2;return;}else{roots=1;root1=root2;}}
}
//#####################################################################
// Function Compute_Roots_Noniterative
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots_Noniterative()
{
    if(c3 == 0){QUADRATIC<T> quadratic(c2,c1,c0);quadratic.Compute_Roots();roots=quadratic.roots;root1=quadratic.root1;root2=quadratic.root2;return;}
    else{
            T one_over_c3=1/c3,a1=c2*one_over_c3,a2=c1*one_over_c3,a3=c0*one_over_c3;
            T Q=(T)one_ninth*a1*a1-(T)one_third*a2,R=a1*((T)one_twenty_seventh*a1*a1-(T)one_sixth*a2)+(T).5*a3;
            T Q_Q_Q=Q*Q*Q,R2_minus_Q3=R*R-Q_Q_Q;       
            if(R2_minus_Q3 <= 0){ // three real roots
                roots=3;
                T theta=acos(R/sqrt(Q_Q_Q)),theta_over_3=(T)one_third*theta,minus_two_sqrt_Q=(T)(-2)*sqrt(Q),minus_a1_over_3=(T)(-one_third)*a1;
                root1=minus_two_sqrt_Q*cos(theta_over_3)+minus_a1_over_3;
                root2=minus_two_sqrt_Q*cos(theta_over_3+(T)two_thirds_pi)+minus_a1_over_3;
                root3=minus_two_sqrt_Q*cos(theta_over_3+(T)four_thirds_pi)+minus_a1_over_3;
                exchange_sort(root1,root2,root3);}      
        else{ // one real root
            roots=1;
            root1=pow(sqrt(R2_minus_Q3)+abs(R),(T)one_third);
            root1+=(T)Q/root1;
            root1*=(R<(T)0)?(T)1:(T)-1;
            root1-=(T)one_third*a1;}}
}
//#####################################################################
// Function Compute_Roots
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots()
{
    root1=0;root2=0;root3=0; // initialize
    if(c3 == 0){QUADRATIC<T> quadratic(c2,c1,c0);quadratic.Compute_Roots();roots=quadratic.roots;root1=quadratic.root1;root2=quadratic.root2;return;}
    else{ // c3 != 0 - cubic - bound the sign changes, i.e. roots
        T bound=(T)1.01*max((T)3*abs(c2/c3),sqrt((T)3*abs(c1/c3)),pow((T)3*abs(c0/c3),(T)1/3));
        Compute_Roots_In_Interval(-bound,bound);}
}
//#####################################################################
// Function Compute_Roots_In_Interval
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Roots_In_Interval(const T xmin,const T xmax)
{
    if(c3 == 0){
        QUADRATIC<T> quadratic(c2,c1,c0);quadratic.Compute_Roots_In_Interval(xmin,xmax);roots=quadratic.roots;root1=quadratic.root1;root2=quadratic.root2;return;}
    else{ // c3 != 0 - cubic
        Compute_Relative_Extrema_Bounding_Sign_Changes_In_Interval(xmin,xmax);
        ITERATIVE_SOLVER<T> iterative_solver;iterative_solver.tolerance=error_tolerance; // use this classes error tolerance
        if(number_of_extrema == 0){
            if((*this)(xmin)*(*this)(xmax) > 0){roots=0;return;}
            else{roots=1;root1=iterative_solver.Bisection_Secant_Root(*this,xmin,xmax);return;}}
        else if(number_of_extrema == 1){roots=2;
            root1=iterative_solver.Bisection_Secant_Root(*this,xmin,extrema1);
            root2=iterative_solver.Bisection_Secant_Root(*this,extrema1,xmax);  
            if(root1 == root2) roots=1;
            return;}
        else{
            roots=3;
            root1=iterative_solver.Bisection_Secant_Root(*this,xmin,extrema1);
            root2=iterative_solver.Bisection_Secant_Root(*this,extrema1,extrema2);   
            root3=iterative_solver.Bisection_Secant_Root(*this,extrema2,xmax); 
            if(root2 == root3) roots=2;
            if(root1 == root2){root2=root3;roots--;}
            return;}}
}
//#####################################################################
// Function Compute_Relative_Extrema
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Relative_Extrema()
{
    QUADRATIC<T> quadratic(3*c3,2*c2,c1);quadratic.Compute_Roots(); // quadratic roots are the extrema of the cubic
    number_of_extrema=quadratic.roots;extrema1=quadratic.root1;extrema2=quadratic.root2;
}
//#####################################################################
// Function Compute_Relative_Extrema_In_Interval
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Relative_Extrema_In_Interval(const T& xmin,const T& xmax)
{
    QUADRATIC<T> quadratic(3*c3,2*c2,c1);quadratic.Compute_Roots_In_Interval(xmin,xmax); // quadratic roots are the extrema of the cubic
    number_of_extrema=quadratic.roots;extrema1=quadratic.root1;extrema2=quadratic.root2;
}
//#####################################################################
// Function Compute_Relative_Extrema_Bounding_Sign_Changes_In_Interval
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Relative_Extrema_Bounding_Sign_Changes_In_Interval(const T& xmin,const T& xmax)
{
    QUADRATIC<T> quadratic(3*c3,2*c2,c1);quadratic.Compute_Roots_In_Interval(xmin,xmax); // quadratic roots are the extrema of the cubic
    number_of_extrema=quadratic.roots;extrema1=quadratic.root1;extrema2=quadratic.root2;
    if(number_of_extrema == 2 && ((*this)(extrema1)*(*this)(extrema2) > 0 || (*this)(extrema2)*(*this)(xmax) > 0)) number_of_extrema=1;
    if(number_of_extrema >= 1 && (*this)(xmin)*(*this)(extrema1) > 0){extrema1=extrema2;number_of_extrema--;}
}
//#####################################################################
// Function Compute_Intervals
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Intervals(const T& a,const T& b,int& intervals,INTERVAL<T>& interval1,INTERVAL<T>& interval2,INTERVAL<T>& interval3)
{
    T ba=b-a,ba2=ba*ba,ba3=ba2*ba;
    T a2=a*a,a3=a2*a;
    CUBIC<T> new_cubic(c3*ba3,(3*c3*a+c2)*ba2,(3*c3*a2+2*c2*a+c1)*ba,c3*a3+c2*a2+c1*a+c0);
    new_cubic.Compute_Intervals(intervals,interval1,interval2,interval3);
    if(intervals){
        interval1=interval1*(b-a)+a;
        if(intervals>1){
            interval2=interval2*(b-a)+a;
            if(intervals>2) interval3=interval3*(b-a)+a;}}
}
//#####################################################################
// Function Compute_Intervals
//#####################################################################
template<class T> void CUBIC<T>::
Compute_Intervals(int& intervals,INTERVAL<T>& interval1,INTERVAL<T>& interval2,INTERVAL<T>& interval3)
{
    // Assume [a,b] are [0,1] and scale appropriately before calling this function
    T g2=c2,g1=c1,g0=c0;
    T h2=3*c3+c2;
    T h1=3*c3+2*c2+c1;
    T h0=c3+c2+c1+c0;
    
    int index=((g0>=0)!=(h0>=0));
    index|=(((g1>=0)!=(h1>=0))<<1);
    index|=(((g2>=0)!=(h2>=0))<<2);
    
    switch(index){
        case 0:
            intervals=0;
            break;
        case 1:
            intervals=1;
            interval1=INTERVAL<T>(0,1);
            break;
        case 2:
        case 6:
            {QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
            if(quadratic.roots!=1){LOG::cout << "2/6 more than one root" << std::endl;return;}
            if((Value(quadratic.root1)>=0)==(g0>=0)) intervals=0;
            else{
                intervals=2;
                interval1=INTERVAL<T>(0,quadratic.root1);
                interval2=INTERVAL<T>(quadratic.root1,1);}
            break;}
        case 3:
        case 7:
            {QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
            if(quadratic.roots!=1){LOG::cout << "3/7 more than one root" << std::endl;return;}
            intervals=1;
            if((Value(quadratic.root1)>=0)==(g0>=0)) interval1=INTERVAL<T>(quadratic.root1,1);
            else interval1=INTERVAL<T>(0,quadratic.root1);
            break;}
        case 4:
            {T inflection_point=-c2/(3*c3);
            QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            if((quadratic(inflection_point)>=0)==(g1>=0)) intervals=0;
            else{
                quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
                if(quadratic.roots!=2){LOG::cout << "4 should be two roots" << std::endl;return;}
                else if((Value(quadratic.root2)>=0)!=(h0>=0)){
                    intervals=2;
                    interval1=INTERVAL<T>(quadratic.root1,quadratic.root2);
                    interval2=INTERVAL<T>(quadratic.root2,1);}
                else if((Value(quadratic.root1)>=0)!=(g0>=0)){
                    intervals=2;
                    interval1=INTERVAL<T>(0,quadratic.root1);
                    interval2=INTERVAL<T>(quadratic.root1,quadratic.root2);}
                else intervals=0;}
            break;}
        case 5:
            {T inflection_point=-c2/(3*c3);
            QUADRATIC<T> quadratic(3*c3,2*c2,c1);
            if((quadratic(inflection_point)>=0)==(g1>=0)){
                intervals=1;
                interval1=INTERVAL<T>(0,1);}
            else{
                quadratic.Compute_Roots_In_Interval(0,1); // quadratic roots are the extrema of the cubic
                T e=Value(quadratic.root1);
                T f=Value(quadratic.root2);
                if(quadratic.roots!=2){
                    intervals=1;
                    interval1=INTERVAL<T>(0,1);}
                else if((e>=0)!=(g0>=0)){
                    if((e>=0)!=(f>=0)){
                        assert((f>=0)!=(h0>=0));
                        intervals=3;
                        interval1=INTERVAL<T>(0,quadratic.root1);
                        interval2=INTERVAL<T>(quadratic.root1,quadratic.root2);
                        interval3=INTERVAL<T>(quadratic.root2,1);}
                    else{
                        intervals=1;
                        interval1=INTERVAL<T>(0,quadratic.root1);}}
                else if((f>=0)!=(h0>=0)){
                    intervals=1;
                    interval1=INTERVAL<T>(quadratic.root2,1);}
                else{
                    assert((e>=0)!=(f>=0));
                    intervals=1;
                    interval1=INTERVAL<T>(quadratic.root1,quadratic.root2);}}
            break;}}
}
//#####################################################################
// Function Coefficients_From_Interpolation
//#####################################################################
template<class T> void CUBIC<T>::
Coefficients_From_Interpolation(T x1,T y1,T x2,T y2,T x3,T y3,T x4,T y4)
{
    T dx1=x1-x4,dx2=x2-x4,dx3=x3-x4,dy1=y1-y4,dy2=y2-y4,dy3=y3-y4;
    T dx12=dx1-dx2,dx13=dx1-dx3,dx23=dx2-dx3,den=1/(dx1*dx2*dx3*dx12*dx13*dx23);
    T sx1=dx1*dx1,sx2=dx2*dx2,sx3=dx3*dx3,n1=dx2*dx3,n2=dx1*dx3,n3=dx1*dx2;
    T k1=n1*dy1,k2=n2*dy2,k3=n3*dy3,m1=dx23*k1,m2=-dx13*k2,m3=dx12*k3;
    T p3=m1+m2+m3;
    T p2=(sx3-sx2)*k1+(sx1-sx3)*k2+(sx2-sx1)*k3;
    T p1=n1*m1+n2*m2+n3*m3;
    c3=den*p3;
    c2=den*(p2-3*p3*x4);
    c1=den*(p1+x4*(3*x4*p3-2*p2));
    c0=y4+den*x4*(x4*(p2-x4*p3)-p1);
}
//#####################################################################
template class CUBIC<float>;
template class CUBIC<double>;
