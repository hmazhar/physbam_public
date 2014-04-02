//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <cmath>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> QUADRATIC<T>::
QUADRATIC()
    :a(1),b(0),c(0),roots(1),root1(0),root2(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> QUADRATIC<T>::
QUADRATIC(const T a_input,const T b_input,const T c_input)
    :a(a_input),b(b_input),c(c_input),roots(0),root1(0),root2(0)
{
}
//#####################################################################
// Function Coefficients_From_Interpolation
//#####################################################################
template<class T> void QUADRATIC<T>::
Coefficients_From_Interpolation(T x1,T y1,T x2,T y2,T x3,T y3)
{
    T dx1=x1-x3,dx2=x2-x3,dy1=y1-y3,dy2=y2-y3,den=1/(dx1*dx2*(dx1-dx2));
    T m1=dy1*dx2,m2=dx1*dy2,p2=m1-m2,p1=dx1*m2-dx2*m1;
    a=den*p2;
    b=den*(p1-2*p2*x3);
    c=den*x3*(x3*p2-p1)+y3;
}
//#####################################################################
// Function Compute_Roots
//#####################################################################
template<class T> void QUADRATIC<T>::
Compute_Roots()
{
    if(a==0){
        if(b==0){
            if(c==0){roots=-1;return;} // function is identically zero - a=b=c=0 - always a root!
            else{roots=0;return;}} // when a=b=0 and c != 0, there are no roots
        else{roots=1;root1=-c/b;return;}} // when a=0 and b != 0, there is one root
    else{ // a != 0
        T d=Discriminant();
        if(d<0){roots=0;return;} // no real roots
        else if(d==0){roots=1;root1=-b/(2*a);return;} // one root
        else{ // d > 0 - two real roots
            T radical;if(b>0) radical=-b-sqrt(d);else radical=-b+sqrt(d);
            roots=2;root1=radical/(2*a);root2=2*c/radical;if(root1>root2) exchange(root1,root2);return;}}
}
//#####################################################################
// Function Compute_Roots_In_Interval
//#####################################################################
template<class T> void QUADRATIC<T>::
Compute_Roots_In_Interval(const T xmin,const T xmax)
{
    Compute_Roots();
    if(roots==1){
        if(root1<xmin || root1>xmax) roots=0;
        else{roots=2;root2=root1;}}
    else if(roots==2){
        if(root2<xmin || root2>xmax) roots--;
        if(root1<xmin || root1>xmax){
            root1=root2;
            roots--;}}
}
template class QUADRATIC<float>;
template class QUADRATIC<double>;
