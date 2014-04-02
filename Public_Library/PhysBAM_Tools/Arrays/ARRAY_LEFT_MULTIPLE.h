//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_LEFT_MULTIPLE
//#####################################################################
#ifndef __ARRAY_LEFT_MULTIPLE__
#define __ARRAY_LEFT_MULTIPLE__

#include <PhysBAM_Tools/Arrays/ARRAY_EXPRESSION.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
namespace PhysBAM{

template<class T1,class T_ARRAY2> class ARRAY_LEFT_MULTIPLE;
template<class T1,class T_ARRAY2> struct IS_ARRAY<ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> > {static const bool value=true;};
template<class T1,class T_ARRAY2> struct IS_ARRAY_VIEW<ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> > {static const bool value=true;};

template<class T1,class T_ARRAY2>
class ARRAY_LEFT_MULTIPLE:public ARRAY_EXPRESSION<typename PRODUCT<T1,typename T_ARRAY2::ELEMENT>::TYPE,ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2>,typename T_ARRAY2::INDEX>
{
    typedef typename T_ARRAY2::ELEMENT T2;
    typedef typename IF<HAS_CHEAP_COPY<T1>::value,const T1,const T1&>::TYPE T1_VIEW; // copy if cheap, otherwise store a reference
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY2>::value,const T_ARRAY2,const T_ARRAY2&>::TYPE T_ARRAY2_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename PRODUCT<T1,T2>::TYPE T_PRODUCT;
public:
    typedef T_PRODUCT ELEMENT;typedef typename T_ARRAY2::INDEX INDEX;

    T1_VIEW c;
    T_ARRAY2_VIEW array;

    ARRAY_LEFT_MULTIPLE(const T1& c,const T_ARRAY2& array)
        :c(c),array(array)
    {}

    INDEX Size() const
    {return array.Size();}

    RANGE<INDEX> Domain_Indices() const
    {return array.Domain_Indices();}

    const T_PRODUCT operator()(const INDEX i) const
    {return c*array(i);}

//#####################################################################
};

template<class T1,class T2,class T_ARRAY2> typename ENABLE_IF<CAN_ASSIGN<typename PRODUCT<T1,typename T_ARRAY2::ELEMENT>::TYPE,typename T_ARRAY2::ELEMENT>::value,ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> >::TYPE
operator*(const T1& c,const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY2::INDEX>& array)
{return ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2>(c,array.Derived());}

template<class T1,class T2,class T_ARRAY2> typename ENABLE_IF<CAN_ASSIGN<typename PRODUCT<T1,typename T_ARRAY2::ELEMENT>::TYPE,typename T_ARRAY2::ELEMENT>::value,ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> >::TYPE
operator/(const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY2::INDEX>& array,const T1& c)
{STATIC_ASSERT(IS_FLOAT_OR_DOUBLE<T1>::value);return ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2>(1/c,array.Derived());}

template<class T1,class T2,class T_ARRAY2> typename ENABLE_IF<CAN_ASSIGN<typename PRODUCT<T1,typename T_ARRAY2::ELEMENT>::TYPE,typename T_ARRAY2::ELEMENT>::value,ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> >::TYPE
operator*(const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY2::INDEX>& array,const T1& c)
{STATIC_ASSERT(IS_SCALAR<T1>::value);return ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2>(c,array.Derived());}

//#####################################################################

template<class T1,class T_ARRAY2> struct PRODUCT<T1,T_ARRAY2,typename IF<IS_SCALAR<T1>::value && IS_ARRAY<T_ARRAY2>::value,void,typename PRODUCT<T1,typename T_ARRAY2::ELEMENT>::TYPE>::TYPE>
{typedef ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> TYPE;};

template<class T1,class T_ARRAY2> struct QUOTIENT<T_ARRAY2,T1,typename IF<IS_FLOAT_OR_DOUBLE<T1>::value && IS_ARRAY<T_ARRAY2>::value,void,typename PRODUCT<T1,typename T_ARRAY2::ELEMENT>::TYPE>::TYPE>
{typedef ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> TYPE;};

template<class T1,class T_ARRAY2> struct PRODUCT<T_ARRAY2,T1,typename IF<IS_SCALAR<T1>::value && IS_ARRAY<T_ARRAY2>::value,void,typename PRODUCT<T1,typename T_ARRAY2::ELEMENT>::TYPE>::TYPE>
{typedef ARRAY_LEFT_MULTIPLE<T1,T_ARRAY2> TYPE;};

//#####################################################################
}
#endif
