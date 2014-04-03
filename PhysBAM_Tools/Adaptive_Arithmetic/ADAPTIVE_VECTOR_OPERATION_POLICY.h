//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_VECTOR_OPERATION_POLICY__
#define __ADAPTIVE_VECTOR_OPERATION_POLICY__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>

namespace PhysBAM{

template<class T1,class T2=ADAPTIVE_DETAIL::NIL,class T3=ADAPTIVE_DETAIL::NIL,class T4=ADAPTIVE_DETAIL::NIL> struct ADAPTIVE_SUM_TYPE
{typedef ADAPTIVE_SUM<ADAPTIVE_SUM<ADAPTIVE_SUM<T1,T2>,T3>,T4> TYPE;};

template<class T1,class T2,class T3> struct ADAPTIVE_SUM_TYPE<T1,T2,T3,ADAPTIVE_DETAIL::NIL>
{typedef ADAPTIVE_SUM<ADAPTIVE_SUM<T1,T2>,T3> TYPE;};

template<class T1,class T2> struct ADAPTIVE_SUM_TYPE<T1,T2,ADAPTIVE_DETAIL::NIL,ADAPTIVE_DETAIL::NIL>
{typedef ADAPTIVE_SUM<T1,T2> TYPE;};

template<class T1> struct ADAPTIVE_SUM_TYPE<T1,ADAPTIVE_DETAIL::NIL,ADAPTIVE_DETAIL::NIL,ADAPTIVE_DETAIL::NIL>
{STATIC_ASSERT(((T1)false));};

template<class EXACT_TYPE,class T,int d>
struct ADAPTIVE_VECTOR_OPERATION_POLICY
{};

template<class T_EXACT,class T>
struct ADAPTIVE_VECTOR_OPERATION_POLICY<T_EXACT,T,2>
{
private:
    typedef typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE T_ADAPTIVE;
public:
    typedef typename ADAPTIVE_SUM_TYPE<T_ADAPTIVE,T_ADAPTIVE>::TYPE SUM;
    typedef ADAPTIVE_DIFFERENCE<ADAPTIVE_PRODUCT<T_ADAPTIVE,T_ADAPTIVE>,ADAPTIVE_PRODUCT<T_ADAPTIVE,T_ADAPTIVE> > DETERMINANT;
};

template<class T_EXACT,class T>
struct ADAPTIVE_VECTOR_OPERATION_POLICY<T_EXACT,T,3>
{
private:
    typedef typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE T_ADAPTIVE;
    typedef typename ADAPTIVE_VECTOR_OPERATION_POLICY<T_EXACT,T,2>::DETERMINANT T_MINOR;
public:
    typedef typename ADAPTIVE_SUM_TYPE<T_ADAPTIVE,T_ADAPTIVE,T_ADAPTIVE>::TYPE SUM;
    typedef typename ADAPTIVE_SUM_TYPE<ADAPTIVE_PRODUCT<T_ADAPTIVE,T_MINOR>,ADAPTIVE_PRODUCT<T_ADAPTIVE,T_MINOR>,ADAPTIVE_PRODUCT<T_ADAPTIVE,T_MINOR> >::TYPE DETERMINANT;
};

template<class T_EXACT,class T>
struct ADAPTIVE_VECTOR_OPERATION_POLICY<T_EXACT,T,4>
{
private:
    typedef typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE T_ADAPTIVE;
public:
    typedef typename ADAPTIVE_SUM_TYPE<T_ADAPTIVE,T_ADAPTIVE,T_ADAPTIVE,T_ADAPTIVE>::TYPE SUM;
};

}
#endif
