//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_PLUS_SCALAR
//#####################################################################
#ifndef __ARRAY_PLUS_SCALAR__
#define __ARRAY_PLUS_SCALAR__

#include <PhysBAM_Tools/Arrays/ARRAY_EXPRESSION.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
namespace PhysBAM{

template<class T1,class T_ARRAY2> class ARRAY_PLUS_SCALAR;
template<class T1,class T_ARRAY2> struct IS_ARRAY<ARRAY_PLUS_SCALAR<T1,T_ARRAY2> > {static const bool value=true;};
template<class T1,class T_ARRAY2> struct IS_ARRAY_VIEW<ARRAY_PLUS_SCALAR<T1,T_ARRAY2> > {static const bool value=true;};

template<class T1,class T_ARRAY2>
class ARRAY_PLUS_SCALAR:public ARRAY_EXPRESSION<typename SUM<T1,typename T_ARRAY2::ELEMENT>::TYPE,ARRAY_PLUS_SCALAR<T1,T_ARRAY2>,typename T_ARRAY2::INDEX>
{
    typedef typename T_ARRAY2::ELEMENT T2;
    typedef typename IF<HAS_CHEAP_COPY<T1>::value,const T1,const T1&>::TYPE T1_VIEW; // copy if cheap, otherwise store a reference
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY2>::value,const T_ARRAY2,const T_ARRAY2&>::TYPE T_ARRAY2_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename SUM<T1,T2>::TYPE T_SUM;
public:
    typedef T_SUM ELEMENT;typedef typename T_ARRAY2::INDEX INDEX;

    T1_VIEW c;
    T_ARRAY2_VIEW array;

    ARRAY_PLUS_SCALAR(const T1& c,const T_ARRAY2& array)
        :c(c),array(array)
    {}

    INDEX Size() const
    {return array.Size();}

    RANGE<INDEX> Domain_Indices() const
    {return array.Domain_Indices();}

    const T_SUM operator()(const INDEX i) const
    {return c+array(i);}

//#####################################################################
};

template<class T1,class T2,class ENABLE=void> struct ARRAY_PLUS_SCALAR_VALID {static const bool value=false;};
template<class T1,class T2> struct ARRAY_PLUS_SCALAR_VALID<T1,T2,typename FIRST<void,typename SUM<T1,T2>::TYPE>::TYPE>
{static const bool value=IS_SAME<T1,T2>::value || IS_SCALAR<T1>::value;};

template<class T1,class T2,class T_ARRAY2> typename ENABLE_IF<ARRAY_PLUS_SCALAR_VALID<T1,T2>::value,ARRAY_PLUS_SCALAR<T1,T_ARRAY2> >::TYPE
operator+(const T1& c,const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY2::INDEX>& array)
{return ARRAY_PLUS_SCALAR<T1,T_ARRAY2>(c,array.Derived());}

template<class T1,class T2,class T_ARRAY2> typename ENABLE_IF<ARRAY_PLUS_SCALAR_VALID<T1,T2>::value,ARRAY_PLUS_SCALAR<T1,T_ARRAY2> >::TYPE
operator+(const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY2::INDEX>& array,const T1& c)
{return ARRAY_PLUS_SCALAR<T1,T_ARRAY2>(c,array.Derived());}

template<class T1,class T2,class T_ARRAY2> typename ENABLE_IF<ARRAY_PLUS_SCALAR_VALID<T1,T2>::value,ARRAY_PLUS_SCALAR<T1,T_ARRAY2> >::TYPE
operator-(const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY2::INDEX>& array,const T1& c)
{return ARRAY_PLUS_SCALAR<T1,T_ARRAY2>(-c,array.Derived());}

//#####################################################################

template<class T1,class T_ARRAY2> struct SUM<T1,T_ARRAY2,typename ENABLE_IF<ARRAY_PLUS_SCALAR_VALID<T1,typename T_ARRAY2::ELEMENT>::value && IS_ARRAY<T_ARRAY2>::value>::TYPE>
{typedef ARRAY_PLUS_SCALAR<T1,T_ARRAY2> TYPE;};

template<class T1,class T_ARRAY2> struct SUM<T_ARRAY2,T1,typename ENABLE_IF<ARRAY_PLUS_SCALAR_VALID<T1,typename T_ARRAY2::ELEMENT>::value && IS_ARRAY<T_ARRAY2>::value>::TYPE>
{typedef ARRAY_PLUS_SCALAR<T1,T_ARRAY2> TYPE;};

template<class T1,class T_ARRAY2> struct DIFFERENCE<T_ARRAY2,T1,typename ENABLE_IF<ARRAY_PLUS_SCALAR_VALID<T1,typename T_ARRAY2::ELEMENT>::value && IS_ARRAY<T_ARRAY2>::value>::TYPE>
{typedef ARRAY_PLUS_SCALAR<T1,T_ARRAY2> TYPE;};

//#####################################################################
}
#endif
