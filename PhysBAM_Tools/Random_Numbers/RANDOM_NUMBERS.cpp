//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Duc Nguyen, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_NUMBERS
//#####################################################################
// See Bratley and Fox. 1988. Algorithm 659: Implementing Sobol's quasirandom sequence generator. ACM Trans. Math. Softw. 14, 88-100.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <limits>
using ::std::log;
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T,class GENERATOR> RANDOM_NUMBERS<T,GENERATOR>::
RANDOM_NUMBERS(const unsigned int seed)
    :gaussian_iset(0)
{
    Set_Seed(seed);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class GENERATOR> RANDOM_NUMBERS<T,GENERATOR>::
~RANDOM_NUMBERS()
{}
//#####################################################################
// Function Set_Seed
//#####################################################################
template<class T,class GENERATOR> void RANDOM_NUMBERS<T,GENERATOR>::
Set_Seed(const unsigned int seed_input)
{
    random_number_generator.Set_Seed(seed_input);
}
//#####################################################################
// Function Get_Uniform_Integer
//#####################################################################
template<class T,class GENERATOR> int RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Integer(const int a,const int b)
{
    return min(b,(int)(a+(b+1-a)*Get_Number())); // in [a,b]
}
//#####################################################################
// Function Get_Uniform_Number
//#####################################################################
template<class T,class GENERATOR> T RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Number(const T a,const T b)
{
    STATIC_ASSERT((NOT<std::numeric_limits<T>::is_integer>::value));
    return a+(b-a)*Get_Number(); // in [a,b)
}
//#####################################################################
// Function Get_Uniform_Vector
//#####################################################################
template<class T,class GENERATOR> template<int d> VECTOR<T,d> RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Vector(const VECTOR<T,d>& v0,const VECTOR<T,d>& v1)
{
    VECTOR<T,d> r;
    for(int i=1;i<=d;i++) r(i)=Get_Uniform_Number(v0(i),v1(i));
    return r;
}
//#####################################################################
// Function Get_Uniform_Vector
//#####################################################################
template<class T,class GENERATOR> template<int d> VECTOR<T,d> RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Vector(const T a,const T b)
{
    VECTOR<T,d> r;
    Fill_Uniform(r,a,b);
    return r;
}
//#####################################################################
// Function Fill_Uniform
//#####################################################################
template<class T,class GENERATOR> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(T& x,const T a,const T b)
{
    x=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform
//#####################################################################
template<class T,class GENERATOR> template<class T_VECTOR> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(VECTOR_BASE<T,T_VECTOR>& v,const T a,const T b)
{
    for(int i=1;i<=v.Size();i++) v(i)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<class T_MATRIX> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(MATRIX_BASE<T,T_MATRIX>& m,const T a,const T b)
{
    for(int i=1;i<=m.Rows();i++) for(int j=1;j<=m.Columns();j++) m(i,j)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<int d> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(DIAGONAL_MATRIX<T,d>& m,const T a,const T b)
{
    for(int i=1;i<=m.Rows();i++) m(i,i)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<int d> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(SYMMETRIC_MATRIX<T,d>& m,const T a,const T b)
{
    for(int i=1;i<=m.Rows();i++) for(int j=1;j<=i;j++) m(i,j)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<int d> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(UPPER_TRIANGULAR_MATRIX<T,d>& m,const T a,const T b)
{
    for(int j=1;j<=d;j++) for(int i=1;i<=j;i++) m(i,j)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform
//#####################################################################
template<class T,class GENERATOR> template<class TV> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(TWIST<TV>& m,const T a,const T b)
{
    Fill_Uniform(m.linear,a,b);
    Fill_Uniform(m.angular,a,b);
}
//#####################################################################
// Function Get_Uniform_Vector
//#####################################################################
template<class T,class GENERATOR> template<class TV> TV RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Vector(const RANGE<TV>& box)
{
    return Get_Uniform_Vector(box.min_corner,box.max_corner);
}
//#####################################################################
// Function Get_Gaussian
//#####################################################################
template<class T,class GENERATOR> T RANDOM_NUMBERS<T,GENERATOR>::
Get_Gaussian()
{
    T fac,rsq,v1,v2;
    if(gaussian_iset==0){
        do{v1=2*Get_Uniform_Number((T)0,(T)1)-1;v2=2*Get_Uniform_Number((T)0,(T)1)-1;rsq=sqr(v1)+sqr(v2);}while(rsq>=1 || rsq==0);
        fac=sqrt(-2*log(rsq)/rsq);gset=v1*fac;gaussian_iset=1;return v2*fac;}
    else{gaussian_iset=0;return gset;}
}
//#####################################################################
// Function Get_Vector_In_Unit_Sphere
//#####################################################################
template<class T,class GENERATOR> template<class TV> TV RANDOM_NUMBERS<T,GENERATOR>::
Get_Vector_In_Unit_Sphere()
{
    for(;;){
        TV v=Get_Uniform_Vector(RANGE<TV>::Centered_Box());
        if(v.Magnitude_Squared()<=1) return v;}
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class T,class GENERATOR> template<class TV> TV RANDOM_NUMBERS<T,GENERATOR>::
Get_Direction()
{
    if(!TV::m) return TV();
    for(;;){
        TV v=Get_Uniform_Vector(RANGE<TV>::Centered_Box());
        typename TV::SCALAR magnitude_squared=v.Magnitude_Squared();
        if(magnitude_squared>0 && magnitude_squared<=1) return v/sqrt(magnitude_squared);}
}
//#####################################################################
// Function Get_Rotation_Helper
//#####################################################################
template<class T>
ROTATION<VECTOR<T,1> > Get_Rotation_Helper(const VECTOR<T,0>&)
{
    return ROTATION<VECTOR<T,1> >();
}
//#####################################################################
// Function Get_Rotation_Helper
//#####################################################################
template<class T>
ROTATION<VECTOR<T,2> > Get_Rotation_Helper(const VECTOR<T,2>& v)
{
    return ROTATION<VECTOR<T,2> >::From_Complex(COMPLEX<T>(v));
}
//#####################################################################
// Function Get_Rotation_Helper
//#####################################################################
template<class T>
ROTATION<VECTOR<T,3> > Get_Rotation_Helper(const VECTOR<T,4>& v)
{
    return ROTATION<VECTOR<T,3> >::From_Quaternion(QUATERNION<T>(v));
}
//#####################################################################
// Function Get_Rotation
//#####################################################################
template<class T,class GENERATOR> template<class TV> ROTATION<TV> RANDOM_NUMBERS<T,GENERATOR>::
Get_Rotation()
{
    return Get_Rotation_Helper(Get_Direction<VECTOR<T,2*TV::m-2> >());
}
//#####################################################################
// Function Get_Frame
//#####################################################################
template<class T,class GENERATOR> template<class TV> FRAME<TV> RANDOM_NUMBERS<T,GENERATOR>::
Get_Frame(const TV& v0,const TV& v1)
{
    TV v=Get_Uniform_Vector(v0,v1);
    return FRAME<TV>(v,Get_Rotation<TV>());
}
//#####################################################################
// Function Get_Twist
//#####################################################################
template<class T,class GENERATOR> template<class TV> TWIST<TV> RANDOM_NUMBERS<T,GENERATOR>::
Get_Twist(const T& a)
{
    TWIST<TV> tw;
    tw.Set_Vector(Get_Uniform_Vector<T,TWIST<TV>::m>(-a,a));
    return tw;
}
//#####################################################################
#define INSTANTIATION_HELPER_V(T,d) \
    template VECTOR<T,d> RANDOM_NUMBERS<T>::Get_Direction<VECTOR<T,d> >(); \
    template VECTOR<T,d> RANDOM_NUMBERS<T>::Get_Uniform_Vector<VECTOR<T,d> >(RANGE<VECTOR<T,d> > const&); \
    template VECTOR<T,d> RANDOM_NUMBERS<T>::Get_Vector_In_Unit_Sphere<VECTOR<T,d> >(); \
    template ROTATION<VECTOR<T,d> > RANDOM_NUMBERS<T>::Get_Rotation<VECTOR<T,d> >(); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<VECTOR<T,d> >(TWIST<VECTOR<T,d> >&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<VECTOR<T,d> >(VECTOR_BASE<T,VECTOR<T,d> >&,T,T)
#define INSTANTIATION_HELPER_V23(T,d) \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<d>(DIAGONAL_MATRIX<T,d>&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<d>(SYMMETRIC_MATRIX<T,d>&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<d>(UPPER_TRIANGULAR_MATRIX<T,d>&,T,T) 
#define IHmn(T,m,n) template void RANDOM_NUMBERS<T>::Fill_Uniform<MATRIX<T,m,n> >(MATRIX_BASE<T,MATRIX<T,m,n> >&,T,T)
#define IHm(T,m) IHmn(T,m,1);IHmn(T,m,2);IHmn(T,m,3);IHmn(T,m,4);IHmn(T,m,5);IHmn(T,m,6)
#define IH(T) IHm(T,1);IHm(T,2);IHm(T,3);IHm(T,4);IHm(T,5);IHm(T,6)
#define INSTANTIATION_HELPER(T) \
    template class RANDOM_NUMBERS<T>; \
    INSTANTIATION_HELPER_V(T,1); \
    INSTANTIATION_HELPER_V(T,2); \
    INSTANTIATION_HELPER_V(T,3); \
    INSTANTIATION_HELPER_V23(T,2); \
    INSTANTIATION_HELPER_V23(T,3); \
    IH(T);                                                              \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<VECTOR<T,0> >(VECTOR_BASE<T,VECTOR<T,0> >&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<VECTOR_ND<T> >(VECTOR_BASE<T,VECTOR_ND<T> >&,T,T); \
    template VECTOR<T,0> RANDOM_NUMBERS<T>::Get_Uniform_Vector<0>(VECTOR<T,0> const&,VECTOR<T,0> const&); \
    template VECTOR<T,1> RANDOM_NUMBERS<T>::Get_Uniform_Vector<1>(VECTOR<T,1> const&,VECTOR<T,1> const&); \
    template VECTOR<T,2> RANDOM_NUMBERS<T>::Get_Uniform_Vector<2>(VECTOR<T,2> const&,VECTOR<T,2> const&); \
    template VECTOR<T,3> RANDOM_NUMBERS<T>::Get_Uniform_Vector<3>(VECTOR<T,3> const&,VECTOR<T,3> const&); \
    template VECTOR<T,4> RANDOM_NUMBERS<T>::Get_Uniform_Vector<4>(VECTOR<T,4> const&,VECTOR<T,4> const&);

INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif
}
