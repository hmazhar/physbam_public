//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_RATIO
//#####################################################################
#ifndef __ARRAY_RATIO__
#define __ARRAY_RATIO__

#include <PhysBAM_Tools/Arrays/ARRAY_EXPRESSION.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
#include <cassert>
namespace PhysBAM{

template<class T_ARRAY1,class T_ARRAY2> class ARRAY_RATIO;
template<class T_ARRAY1,class T_ARRAY2> struct IS_ARRAY<ARRAY_RATIO<T_ARRAY1,T_ARRAY2> > {static const bool value=true;};
template<class T_ARRAY1,class T_ARRAY2> struct IS_ARRAY_VIEW<ARRAY_RATIO<T_ARRAY1,T_ARRAY2> > {static const bool value=true;};

template<class T_ARRAY1,class T_ARRAY2>
class ARRAY_RATIO:public ARRAY_EXPRESSION<typename T_ARRAY1::ELEMENT,ARRAY_RATIO<T_ARRAY1,T_ARRAY2>,typename T_ARRAY1::INDEX>
{
    typedef typename T_ARRAY1::ELEMENT T1;typedef typename T_ARRAY2::ELEMENT T2;
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY1>::value,const T_ARRAY1,const T_ARRAY1&>::TYPE T_ARRAY1_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY2>::value,const T_ARRAY2,const T_ARRAY2&>::TYPE T_ARRAY2_VIEW;
    typedef typename QUOTIENT<T1,T2>::TYPE T_RATIO;
public:
    typedef T_RATIO ELEMENT;typedef typename T_ARRAY1::INDEX INDEX;

    T_ARRAY1_VIEW array1;
    T_ARRAY2_VIEW array2;

    ARRAY_RATIO(const T_ARRAY1& array1,const T_ARRAY2& array2)
        :array1(array1),array2(array2)
    {}

    INDEX Size() const
    {INDEX size=array1.Size();assert(size==array2.Size());return size;}

    RANGE<INDEX> Domain_Indices() const
    {RANGE<INDEX> domain_indices=array1.Domain_Indices();assert(domain_indices==array2.Domain_Indices());return domain_indices;}

    const T_RATIO operator()(const INDEX i) const
    {return array1(i)/array2(i);}

//#####################################################################
};

template<class T1,class T2,class T_ARRAY1,class T_ARRAY2> ARRAY_RATIO<T_ARRAY1,T_ARRAY2>
operator/(const ARRAY_BASE<T1,T_ARRAY1,typename T_ARRAY1::INDEX>& array1,const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY1::INDEX>& array2)
{return ARRAY_RATIO<T_ARRAY1,T_ARRAY2>(array1.Derived(),array2.Derived());}

//#####################################################################

template<class T_ARRAY1,class T_ARRAY2> struct QUOTIENT<T_ARRAY1,T_ARRAY2,typename ENABLE_IF<IS_ARRAY<T_ARRAY1>::value && IS_ARRAY<T_ARRAY2>::value>::TYPE>
{typedef ARRAY_RATIO<T_ARRAY1,T_ARRAY2> TYPE;};

//#####################################################################

}
#endif
