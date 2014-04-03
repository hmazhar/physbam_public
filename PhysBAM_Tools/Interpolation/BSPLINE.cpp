//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Interpolation/BSPLINE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> BSPLINE<T,T2>::
BSPLINE(const ARRAY<T>& control_points_times,const ARRAY<T2>& control_points,const int order)
    :control_points_times(control_points_times),control_points(control_points),k(order+1),closed(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> BSPLINE<T,T2>::
~BSPLINE()
{
}
//#####################################################################
// Function Evaluate
//#####################################################################
template<class T,class T2> T2 BSPLINE<T,T2>::
Evaluate(const T t)
{
    assert(Start_Time()<=t && t<End_Time());
    T2 sum=T2();
    for(int i=1;i<=control_points.m-k;i++){sum+=Basis_Function(i,k,t)*control_points(i+k/2);}
    return sum;
}
//#####################################################################
// Function Clamped_Evaluate
//#####################################################################
template<class T,class T2> T2 BSPLINE<T,T2>::
Clamped_Evaluate(const T t)
{
    T2 sum=T2();
    T clamped_t=clamp(t,Start_Time(),End_Time()-(T)1e-7);
    for(int i=1;i<=control_points.m-k;i++){sum+=Basis_Function(i,k,clamped_t)*control_points(i+k/2);}
    return sum;
}
//#####################################################################
// Function Basis_Function
//#####################################################################
template<class T,class T2> T BSPLINE<T,T2>::
Basis_Function(const int i,const int k,const T t)
{
    T t_i=control_points_times(i),t_i_plus_1=control_points_times(i+1);
    if(k==1){if(t_i<=t && t<t_i_plus_1) return 1;else return 0;}
    T t_i_plus_k_minus_1=control_points_times(i+k-1),t_i_plus_k=control_points_times(i+k);
    T result=0;
    if(t_i_plus_k_minus_1!=t_i) result+=((t-t_i)/(t_i_plus_k_minus_1-t_i))*Basis_Function(i,k-1,t);
    if(t_i_plus_k!=t_i_plus_1) result+=((t_i_plus_k-t)/(t_i_plus_k-t_i_plus_1))*Basis_Function(i+1,k-1,t);
    return result;
}
//#####################################################################
// Function Clamp_End_Points
//#####################################################################
template<class T,class T2> void BSPLINE<T,T2>::
Clamp_End_Points()
{
    assert(!closed);assert(k>=2);
    int number_of_knots=control_points_times.m;int new_points=k-1;
    ARRAY<T> new_control_points_times;new_control_points_times.Resize(number_of_knots+2*new_points);
    ARRAY<T2> new_control_points;new_control_points.Resize(number_of_knots+2*new_points);
    for(int i=1;i<=number_of_knots;i++){new_control_points_times(i+new_points)=control_points_times(i);new_control_points(i+new_points)=control_points(i);}
    for(int i=1;i<=new_points;i++){
        new_control_points_times(i)=control_points_times(1);new_control_points_times(number_of_knots+new_points+i)=control_points_times(number_of_knots);
        new_control_points(i)=control_points(1);new_control_points(number_of_knots+new_points+i)=control_points(number_of_knots);}
    control_points_times=new_control_points_times;control_points=new_control_points;
}
//#####################################################################
// Function Create_Closed_Points
//#####################################################################
template<class T,class T2> void BSPLINE<T,T2>::
Create_Closed_Points()
{
    control_points.Resize(control_points.m+2*(k-1));control_points_times.Resize(control_points_times.m+2*(k-1));
    for(int i=1;i<=2*(k-1);i++){
        control_points(control_points.m-2*(k-1)+i)=control_points(i+1);
        control_points_times(control_points_times.m-2*(k-1)+i)=(control_points_times(i+1)-control_points_times(i)+control_points_times(control_points_times.m-2*(k-1)+i-1));}
    closed=true;
}
//#####################################################################
// Function Normalize_Control_Points
//#####################################################################
template<class T,class T2> void BSPLINE<T,T2>::
Normalize_Control_Points()
{
    PHYSBAM_FATAL_ERROR(); // TODO: This form of sort does not exist.
//     ARRAY<T>::sort(control_points_times,control_points);
//     T offset=Start_Time(),total=Range();
//     for(int i=1;i<=control_points_times.m;i++) control_points_times(i)=(control_points_times(i)-offset)/total;
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
//#####################################################################
// Function Print_Control_Points_And_Times
//#####################################################################
template<class T,class T2> void BSPLINE<T,T2>::
Print_Control_Points_And_Times()
{
    for(int i=1;i<=control_points.m;i++) LOG::cout<<"Time: "<<control_points_times(i)<<", Control Point: "<<control_points(i)<<std::endl;
}
#endif
template class BSPLINE<float,VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BSPLINE<double,VECTOR<double,3> >;
#endif
