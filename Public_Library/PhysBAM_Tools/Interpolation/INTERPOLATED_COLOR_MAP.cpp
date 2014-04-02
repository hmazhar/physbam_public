//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> INTERPOLATED_COLOR_MAP<T>::
INTERPOLATED_COLOR_MAP()
    :mn(0),mx(1),use_log(false),colors(*new INTERPOLATION_CURVE<T,VECTOR<T,3> >)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> INTERPOLATED_COLOR_MAP<T>::
~INTERPOLATED_COLOR_MAP()
{
    delete &colors;
}
//#####################################################################
// Function Initialize_Colors
//#####################################################################
template<class T> void INTERPOLATED_COLOR_MAP<T>::
Initialize_Colors(T min_value,T max_value,bool log_scale,bool reverse,bool two_sets)
{
    mn=min_value;
    mx=max_value;
    use_log=log_scale;
    if(use_log){
        mn=std::log(std::abs(mn));
        mx=std::log(std::abs(mx));}
    PHYSBAM_ASSERT(mn<mx);
    T a=mn,d=(mx-mn)/(6*(1+two_sets));
    if(reverse){a=mx,d=-d;}
    colors.Add_Control_Point(a+0*d,VECTOR<T,3>(1,1,1));
    colors.Add_Control_Point(a+1*d,VECTOR<T,3>(1,0,0));
    colors.Add_Control_Point(a+2*d,VECTOR<T,3>(1,1,0));
    colors.Add_Control_Point(a+3*d,VECTOR<T,3>(0,1,0));
    colors.Add_Control_Point(a+4*d,VECTOR<T,3>(0,1,1));
    colors.Add_Control_Point(a+5*d,VECTOR<T,3>(0,0,1));
    if(two_sets){
        colors.Add_Control_Point(a+6*d,VECTOR<T,3>(.3,0,.3));
        colors.Add_Control_Point(a+7*d,VECTOR<T,3>(.3,0,0));
        colors.Add_Control_Point(a+8*d,VECTOR<T,3>(.3,.3,0));
        colors.Add_Control_Point(a+9*d,VECTOR<T,3>(0,.3,0));
        colors.Add_Control_Point(a+10*d,VECTOR<T,3>(0,.3,.3));
        colors.Add_Control_Point(a+11*d,VECTOR<T,3>(0,0,.3));}
    colors.Add_Control_Point(a+6*(1+two_sets)*d,VECTOR<T,3>());
}
//#####################################################################
// Function operator ()
//#####################################################################
template<class T> VECTOR<T,3> INTERPOLATED_COLOR_MAP<T>::
operator()(T x) const
{
    if(use_log){if(x==0) x=mn;else x=std::log(std::abs(x));}
    return colors.Value(x);
}
template class INTERPOLATED_COLOR_MAP<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERPOLATED_COLOR_MAP<double>;
#endif
