//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_SIGNED_VOLUME__
#define __ADAPTIVE_SIGNED_VOLUME__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_DETERMINANT.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_ADAPTIVE,int d> struct IS_ADAPTIVE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,d> > {static const bool value=true;};

template<class T_ADAPTIVE>
class ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,1>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,1>,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,1> >
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,1> >::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a11,a21;
    ADAPTIVE_SIGNED_VOLUME();
    ADAPTIVE_SIGNED_VOLUME& operator=(const ADAPTIVE_SIGNED_VOLUME&);

public:
    ADAPTIVE_SIGNED_VOLUME(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a21_input)
        :a11(a11_input),a21(a21_input)
    {}

public:
    WRAPPED_TYPE Wrapped_Implementation() const
    {return Adaptive_Determinant<void>(a11-a21);}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,2>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,2>,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,2> >
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,2> >::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a11,a12,a21,a22,a31,a32;
    ADAPTIVE_SIGNED_VOLUME();
    ADAPTIVE_SIGNED_VOLUME& operator=(const ADAPTIVE_SIGNED_VOLUME&);

public:
    ADAPTIVE_SIGNED_VOLUME(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input,const T_ADAPTIVE& a31_input,const T_ADAPTIVE& a32_input)
        :a11(a11_input),a12(a12_input),a21(a21_input),a22(a22_input),a31(a31_input),a32(a32_input)
    {}

public:
    WRAPPED_TYPE Wrapped_Implementation() const
    {return Adaptive_Determinant<void>(a11-a31,a12-a32,a21-a31,a22-a32);}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,3>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,3>,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,3> >
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,3> >::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a11,a12,a13,a21,a22,a23,a31,a32,a33,a41,a42,a43;
    ADAPTIVE_SIGNED_VOLUME();
    ADAPTIVE_SIGNED_VOLUME& operator=(const ADAPTIVE_SIGNED_VOLUME&);

public:
    ADAPTIVE_SIGNED_VOLUME(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a13_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input,const T_ADAPTIVE& a23_input,
        const T_ADAPTIVE& a31_input,const T_ADAPTIVE& a32_input,const T_ADAPTIVE& a33_input,const T_ADAPTIVE& a41_input,const T_ADAPTIVE& a42_input,const T_ADAPTIVE& a43_input)
        :a11(a11_input),a12(a12_input),a13(a13_input),a21(a21_input),a22(a22_input),a23(a23_input),a31(a31_input),a32(a32_input),a33(a33_input),a41(a41_input),a42(a42_input),a43(a43_input)
    {}

public:
    WRAPPED_TYPE Wrapped_Implementation() const
    {return Adaptive_Determinant<void>(a11-a41,a12-a42,a13-a43,a21-a41,a22-a42,a23-a43,a31-a41,a32-a42,a33-a43);}
//#####################################################################
};

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const T& t11,const T& t21)
{return ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>(t11,t21);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const VECTOR<T,1>& t1,const VECTOR<T,1>& t2)
{return Adaptive_Signed_Volume<T_EXACT>(t1[1],t2[1]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const VECTOR<T,1>& t1,const VECTOR<VECTOR<T,1>,1>& t2)
{return Adaptive_Signed_Volume<T_EXACT>(t1,t2[1]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const VECTOR<VECTOR<T,1>,2>& t)
{return Adaptive_Signed_Volume<T_EXACT,T>(t[1],t[2]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const T& t11,const T& t12,const T& t21,const T& t22,const T& t31,const T& t32)
{return ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>(t11,t12,t21,t22,t31,t32);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const VECTOR<T,2>& t1,const VECTOR<T,2>& t2,const VECTOR<T,2>& t3)
{return Adaptive_Signed_Volume<T_EXACT>(t1[1],t1[2],t2[1],t2[2],t3[1],t3[2]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const VECTOR<T,2>& t1,const VECTOR<VECTOR<T,2>,2>& t23)
{return Adaptive_Signed_Volume<T_EXACT>(t1,t23[1],t23[2]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const VECTOR<VECTOR<T,2>,3>& t)
{return Adaptive_Signed_Volume<T_EXACT>(t[1],t[2],t[3]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const T& t11,const T& t12,const T& t13,const T& t21,const T& t22,const T& t23,const T& t31,const T& t32,const T& t33,const T& t41,const T& t42,const T& t43)
{return ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>(t11,t12,t13,t21,t22,t23,t31,t32,t33,t41,t42,t43);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const VECTOR<T,3>& t1,const VECTOR<T,3>& t2,const VECTOR<T,3>& t3,const VECTOR<T,3>& t4)
{return Adaptive_Signed_Volume<T_EXACT>(t1[1],t1[2],t1[3],t2[1],t2[2],t2[3],t3[1],t3[2],t3[3],t4[1],t4[2],t4[3]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const VECTOR<T,3>& t1,const VECTOR<VECTOR<T,3>,3>& t234)
{return Adaptive_Signed_Volume<T_EXACT>(t1,t234[1],t234[2],t234[3]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const VECTOR<VECTOR<T,3>,4>& t)
{return Adaptive_Signed_Volume<T_EXACT>(t[1],t[2],t[3],t[4]);}

}
using ADAPTIVE_DETAIL::ADAPTIVE_SIGNED_VOLUME;
using ADAPTIVE_DETAIL::Adaptive_Signed_Volume;
}
#endif
