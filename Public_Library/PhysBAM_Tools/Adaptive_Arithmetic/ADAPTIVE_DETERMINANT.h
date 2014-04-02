//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_DETERMINANT__
#define __ADAPTIVE_DETERMINANT__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OP.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_ADAPTIVE,int d> struct IS_ADAPTIVE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,d> > {static const bool value=true;};

template<class T_ADAPTIVE>
class ADAPTIVE_DETERMINANT<T_ADAPTIVE,1>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,1>,T_ADAPTIVE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT,T_ADAPTIVE>::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a11;
    ADAPTIVE_DETERMINANT();
    void operator=(const ADAPTIVE_DETERMINANT&);

public:
    ADAPTIVE_DETERMINANT(const T_ADAPTIVE& a11_input)
        :a11(a11_input)
    {}

    WRAPPED_TYPE Wrapped_Implementation() const
    {return a11;}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_DETERMINANT<T_ADAPTIVE,2>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,2>,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,2>::DETERMINANT>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,2>::DETERMINANT>::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a11,a12,a21,a22;
    ADAPTIVE_DETERMINANT();
    void operator=(const ADAPTIVE_DETERMINANT&);

public:
    ADAPTIVE_DETERMINANT(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input)
        :a11(a11_input),a12(a12_input),a21(a21_input),a22(a22_input)
    {}

    WRAPPED_TYPE Wrapped_Implementation() const
    {return (a11*a22-a12*a21);}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_DETERMINANT<T_ADAPTIVE,3>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,3>,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,3>::DETERMINANT>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,3>::DETERMINANT>::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a11,a12,a13,a21,a22,a23,a31,a32,a33;
    ADAPTIVE_DETERMINANT();
    void operator=(const ADAPTIVE_DETERMINANT&);

public:
    ADAPTIVE_DETERMINANT(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a13_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input,
        const T_ADAPTIVE& a23_input,const T_ADAPTIVE& a31_input,const T_ADAPTIVE& a32_input,const T_ADAPTIVE& a33_input)
        :a11(a11_input),a12(a12_input),a13(a13_input),a21(a21_input),a22(a22_input),a23(a23_input),a31(a31_input),a32(a32_input),a33(a33_input)
    {}

    WRAPPED_TYPE Wrapped_Implementation() const
    {return (a11*(a22*a33-a23*a32)+a12*(a23*a31-a21*a33)+a13*(a21*a32-a22*a31));}
//#####################################################################
};

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Determinant(const T& t11)
{return ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>(t11);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Determinant(const VECTOR<T,1>& t1)
{return Adaptive_Determinant<T_EXACT>(t1[1]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Determinant(const VECTOR<VECTOR<T,1>,1>& t)
{return Adaptive_Determinant<T_EXACT>(t[1]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Determinant(const T& t11,const T& t12,const T& t21,const T& t22)
{return ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>(t11,t12,t21,t22);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Determinant(const VECTOR<T,2>& t1,const VECTOR<T,2>& t2)
{return Adaptive_Determinant<T_EXACT>(t1[1],t2[1],t1[2],t2[2]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Determinant(const VECTOR<VECTOR<T,2>,2>& t)
{return Adaptive_Determinant<T_EXACT>(t[1],t[2]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Determinant(const T& t11,const T& t12,const T& t13,const T& t21,const T& t22,const T& t23,const T& t31,const T& t32,const T& t33)
{return ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>(t11,t12,t13,t21,t22,t23,t31,t32,t33);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Determinant(const VECTOR<T,3>& t1,const VECTOR<T,3>& t2,const VECTOR<T,3>& t3)
{return Adaptive_Determinant<T_EXACT>(t1[1],t2[1],t3[1],t1[2],t2[2],t3[2],t1[3],t2[3],t3[3]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Determinant(const VECTOR<VECTOR<T,3>,3>& t)
{return Adaptive_Determinant<T_EXACT>(t[1],t[2],t[3]);}

}
using ADAPTIVE_DETAIL::ADAPTIVE_DETERMINANT;
using ADAPTIVE_DETAIL::Adaptive_Determinant;
}
#endif
