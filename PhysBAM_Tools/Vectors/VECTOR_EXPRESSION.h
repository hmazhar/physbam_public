//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_EXPRESSION
//#####################################################################
#ifndef __VECTOR_EXPRESSION__
#define __VECTOR_EXPRESSION__

#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_BASE.h>
#include <cassert>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

template<class T,class T_EXPRESSION>
class VECTOR_EXPRESSION:public VECTOR_BASE<T,T_EXPRESSION>
{
};

template<class T_VECTOR1,class T_VECTOR2>
class VECTOR_SUM:public VECTOR_EXPRESSION<typename SUM<typename T_VECTOR1::SCALAR,typename T_VECTOR2::SCALAR>::TYPE,VECTOR_SUM<T_VECTOR1,T_VECTOR2> >
{
public:
    typedef VECTOR_EXPRESSION<typename SUM<typename T_VECTOR1::SCALAR,typename T_VECTOR2::SCALAR>::TYPE,VECTOR_SUM<T_VECTOR1,T_VECTOR2> > BASE;
    typedef typename T_VECTOR1::SCALAR T1;typedef typename T_VECTOR2::SCALAR T2;
    typedef typename SUM<T1,T2>::TYPE SCALAR;
    typedef SCALAR ELEMENT;

    const T_VECTOR1& array1;
    const T_VECTOR2& array2;

    VECTOR_SUM(const T_VECTOR1& array1,const T_VECTOR2& array2)
        :array1(array1),array2(array2)
    {}

    VECTOR_SUM(const VECTOR_SUM& v)
        :BASE(),array1(v.array1),array2(v.array2)
    {}

    int Size() const
    {assert(array1.Size()==array2.Size());return array1.Size();}

    const SCALAR operator()(const int i) const
    {return array1(i)+array2(i);}

//#####################################################################
};

template<class T_VECTOR1,class T_VECTOR2>
class VECTOR_DIFFERENCE:public VECTOR_EXPRESSION<typename SUM<typename T_VECTOR1::SCALAR,typename T_VECTOR2::SCALAR>::TYPE,VECTOR_DIFFERENCE<T_VECTOR1,T_VECTOR2> >
{
public:
    typedef VECTOR_EXPRESSION<typename SUM<typename T_VECTOR1::SCALAR,typename T_VECTOR2::SCALAR>::TYPE,VECTOR_DIFFERENCE<T_VECTOR1,T_VECTOR2> > BASE;
    typedef typename T_VECTOR1::SCALAR T1;typedef typename T_VECTOR2::SCALAR T2;
    typedef typename DIFFERENCE<T1,T2>::TYPE SCALAR;
    typedef SCALAR ELEMENT;

    const T_VECTOR1& array1;
    const T_VECTOR2& array2;

    VECTOR_DIFFERENCE(const T_VECTOR1& array1,const T_VECTOR2& array2)
        :array1(array1),array2(array2)
    {}

    VECTOR_DIFFERENCE(const VECTOR_DIFFERENCE& v)
        :BASE(),array1(v.array1),array2(v.array2)
    {}

    int Size() const
    {assert(array1.Size()==array2.Size());return array1.Size();}

    const SCALAR operator()(const int i) const
    {return array1(i)-array2(i);}

//#####################################################################
};

template<class T1,class T_VECTOR>
class VECTOR_SCALE:public VECTOR_EXPRESSION<typename SUM<T1,typename T_VECTOR::SCALAR>::TYPE,VECTOR_SCALE<T1,T_VECTOR> >
{
public:
    typedef VECTOR_EXPRESSION<typename SUM<T1,typename T_VECTOR::SCALAR>::TYPE,VECTOR_SCALE<T1,T_VECTOR> > BASE;
    typedef typename T_VECTOR::SCALAR T2;
    typedef typename PRODUCT<T1,T2>::TYPE SCALAR;
    typedef SCALAR ELEMENT;

    const T_VECTOR& array;
    T1 scalar;

    VECTOR_SCALE(const T1& scalar,const T_VECTOR& array)
        :array(array),scalar(scalar)
    {}

    VECTOR_SCALE(const VECTOR_SCALE& v)
        :BASE(),array(v.array),scalar(v.scalar)
    {}

    int Size() const
    {return array.Size();}

    const SCALAR operator()(const int i) const
    {return scalar*array(i);}

//#####################################################################
};

template<class T_VECTOR>
class VECTOR_NEGATION:public VECTOR_EXPRESSION<typename T_VECTOR::SCALAR,VECTOR_NEGATION<T_VECTOR> >
{
public:
    typedef VECTOR_EXPRESSION<typename T_VECTOR::SCALAR,VECTOR_NEGATION<T_VECTOR> > BASE;
    typedef typename NEGATION<typename T_VECTOR::SCALAR>::TYPE SCALAR;
    typedef SCALAR ELEMENT;

    const T_VECTOR& array;

    VECTOR_NEGATION(const T_VECTOR& array)
        :array(array)
    {}

    VECTOR_NEGATION(const VECTOR_NEGATION& v)
        :BASE(),array(v.array)
    {}

    int Size() const
    {return array.Size();}

    const SCALAR operator()(const int i) const
    {return -array(i);}

//#####################################################################
};

template<class T1,class T2,class T_VECTOR1,class T_VECTOR2> VECTOR_SUM<T_VECTOR1,T_VECTOR2>
operator+(const VECTOR_BASE<T1,T_VECTOR1>& array1,const VECTOR_BASE<T2,T_VECTOR2>& array2)
{array1.Static_Assert_Not_Small();array2.Static_Assert_Not_Small();return VECTOR_SUM<T_VECTOR1,T_VECTOR2>(array1.Derived(),array2.Derived());}

template<class T1,class T2,class T_VECTOR1,class T_VECTOR2> VECTOR_DIFFERENCE<T_VECTOR1,T_VECTOR2>
operator-(const VECTOR_BASE<T1,T_VECTOR1>& array1,const VECTOR_BASE<T2,T_VECTOR2>& array2)
{array1.Static_Assert_Not_Small();array2.Static_Assert_Not_Small();return VECTOR_DIFFERENCE<T_VECTOR1,T_VECTOR2>(array1.Derived(),array2.Derived());}

template<class T,class T_VECTOR> VECTOR_NEGATION<T_VECTOR>
operator-(const VECTOR_BASE<T,T_VECTOR>& array)
{array.Static_Assert_Not_Small();return VECTOR_NEGATION<T_VECTOR>(array.Derived());}

template<class T,class T_VECTOR2> VECTOR_SCALE<T,T_VECTOR2>
operator*(const T& c,const VECTOR_BASE<T,T_VECTOR2>& array)
{array.Static_Assert_Not_Small();return VECTOR_SCALE<T,T_VECTOR2>(c,array.Derived());}

template<class T,class T_VECTOR2> VECTOR_SCALE<T,T_VECTOR2>
operator/(const VECTOR_BASE<T,T_VECTOR2>& array,const T& c)
{array.Static_Assert_Not_Small();return VECTOR_SCALE<T,T_VECTOR2>(1/c,array.Derived());}

template<class T,class T_VECTOR2> VECTOR_SCALE<T,T_VECTOR2>
operator*(const VECTOR_BASE<T,T_VECTOR2>& array,const T& c)
{array.Static_Assert_Not_Small();return VECTOR_SCALE<T,T_VECTOR2>(c,array.Derived());}

//#####################################################################

template<class T_VECTOR1,class T_VECTOR2> struct SUM<T_VECTOR1,T_VECTOR2,typename ENABLE_IF<INEFFICIENT_VECTOR<T_VECTOR1>::value && INEFFICIENT_VECTOR<T_VECTOR2>::value>::TYPE>
{typedef VECTOR_SUM<T_VECTOR1,T_VECTOR2> TYPE;};

template<class T_VECTOR1,class T_VECTOR2> struct DIFFERENCE<T_VECTOR1,T_VECTOR2,typename ENABLE_IF<INEFFICIENT_VECTOR<T_VECTOR1>::value && INEFFICIENT_VECTOR<T_VECTOR2>::value>::TYPE>
{typedef VECTOR_DIFFERENCE<T_VECTOR1,T_VECTOR2> TYPE;};

template<class T_VECTOR> struct NEGATION<T_VECTOR,typename ENABLE_IF<INEFFICIENT_VECTOR<T_VECTOR>::value>::TYPE>
{typedef VECTOR_NEGATION<T_VECTOR> TYPE;};

template<class T1,class T_VECTOR2> struct PRODUCT<T1,T_VECTOR2,typename ENABLE_IF<INEFFICIENT_VECTOR<T_VECTOR2>::value && IS_SCALAR<T1>::value>::TYPE>
{typedef VECTOR_SCALE<T1,T_VECTOR2> TYPE;};

template<class T1,class T_VECTOR2> struct QUOTIENT<T_VECTOR2,T1,typename ENABLE_IF<INEFFICIENT_VECTOR<T_VECTOR2>::value && IS_SCALAR<T1>::value>::TYPE>
{typedef VECTOR_SCALE<T1,T_VECTOR2> TYPE;};

template<class T1,class T_VECTOR2> struct PRODUCT<T_VECTOR2,T1,typename ENABLE_IF<INEFFICIENT_VECTOR<T_VECTOR2>::value && IS_SCALAR<T1>::value>::TYPE>
{typedef VECTOR_SCALE<T1,T_VECTOR2> TYPE;};

//#####################################################################

}
#endif
