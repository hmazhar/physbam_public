//#####################################################################
// Copyright 2005, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Interpolation/BSPLINE_QUATERNION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> BSPLINE_QUATERNION<T>::
BSPLINE_QUATERNION(const ARRAY<T>& control_points_times,const ARRAY<ROTATION<TV> >& control_points,const int order)
    :BSPLINE<T,ROTATION<TV> >(control_points_times,control_points,order)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> BSPLINE_QUATERNION<T>::
~BSPLINE_QUATERNION()
{
}
//#####################################################################
// Function Evaluate
//#####################################################################
template<class T> ROTATION<VECTOR<T,3> > BSPLINE_QUATERNION<T>::
Evaluate(const T t)
{
    ROTATION<TV>  total;
    int l=1;
    for(int i=1;i<=control_points.m-k;i++) if(control_points_times(i)>t){l=i-k;break;}
    if(l<1){l=1;total=ROTATION<TV>::From_Rotation_Vector(control_points(1).Rotation_Vector()*Quaternion_Basis_Function(1,k,t));assert(!closed);}else{total=control_points(l);}
    for(int i=l+1;i<=control_points.m-k && control_points_times(i)<t;i++) total*=ROTATION<TV>::From_Rotation_Vector(Omega(i).Rotation_Vector()*Quaternion_Basis_Function(i,k,t));
    return total;
}
//#####################################################################
// Function Quaternion_Basis_Function
//#####################################################################
template<class T> T BSPLINE_QUATERNION<T>::
Quaternion_Basis_Function(const int i,const int k,const T t)
{
    T sum=0;
    for(int j=i;j<=control_points.m-k;j++) sum+=BSPLINE<T,ROTATION<TV> >::Basis_Function(j,k,t);
    return sum;
}
//#####################################################################
// Function Omega
//#####################################################################
template<class T> ROTATION<VECTOR<T,3> > BSPLINE_QUATERNION<T>::
Omega(const int i)
{
    return control_points(i-1).Inverse()*control_points(i);
}
//#####################################################################
// Function Create_Closed_Points
//#####################################################################
template<class T> void BSPLINE_QUATERNION<T>::
Create_Closed_Points()
{
    BSPLINE<T,ROTATION<TV> >::Create_Closed_Points();
    for(int i=k-1;i>=1;i--){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Comparing "<<i<<": "<<control_points(i)<<" and "<<i+1<<": "<<control_points(i+1)<<std::endl;
#endif
        control_points(i)=ROTATION<TV>::Switch_Hemisphere(control_points(i+1),control_points(i));
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"   result: "<<control_points(i)<<std::endl;
#endif
    }
    if(k>2) for(int i=control_points.m-(k-2);i<control_points.m;i++){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Comparing "<<i<<": "<<control_points(i)<<" and "<<i+1<<": "<<control_points(i+1)<<std::endl;
#endif
        control_points(i)=ROTATION<TV>::Switch_Hemisphere(control_points(i+1),control_points(i));
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"   result: "<<control_points(i)<<std::endl;
#endif
    }
}
//#####################################################################
// Function Quaternion_Check
//#####################################################################
template<class T> void BSPLINE_QUATERNION<T>::
Quaternion_Check()
{
    for(int i=2;i<=control_points.m;i++) control_points(i)=ROTATION<TV>::Switch_Hemisphere(control_points(i-1),control_points(i));
}
