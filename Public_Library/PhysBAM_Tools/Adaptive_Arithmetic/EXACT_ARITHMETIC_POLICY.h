//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXACT_ARITHMETIC_POLICY
//#####################################################################
#ifndef __EXACT_ARITHMETIC_POLICY__
#define __EXACT_ARITHMETIC_POLICY__

namespace PhysBAM{

template<class T> class EXACT_FLOAT;
template<class T> class EXACT_RATIONAL;

template<class T1, class T2> struct EXACT_ARITHMETIC_POLICY;

//#####################################################################
// EXACT_FLOAT specializations
//#####################################################################
template<class T>
struct EXACT_ARITHMETIC_POLICY<EXACT_FLOAT<T>,T>
{
    typedef EXACT_FLOAT<T> SUM_TYPE;
    typedef EXACT_FLOAT<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
template<class T>
struct EXACT_ARITHMETIC_POLICY<T,EXACT_FLOAT<T> >
{
    typedef EXACT_FLOAT<T> SUM_TYPE;
    typedef EXACT_FLOAT<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
template<class T>
struct EXACT_ARITHMETIC_POLICY<EXACT_FLOAT<T>,EXACT_FLOAT<T> >
{
    typedef EXACT_FLOAT<T> SUM_TYPE;
    typedef EXACT_FLOAT<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
//#####################################################################
// EXACT_RATIONAL specializations
//#####################################################################
template<class T>
struct EXACT_ARITHMETIC_POLICY<EXACT_RATIONAL<T>,T>
{
    typedef EXACT_RATIONAL<T> SUM_TYPE;
    typedef EXACT_RATIONAL<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
template<class T>
struct EXACT_ARITHMETIC_POLICY<EXACT_RATIONAL<T>,EXACT_FLOAT<T> >
{
    typedef EXACT_RATIONAL<T> SUM_TYPE;
    typedef EXACT_RATIONAL<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
template<class T>
struct EXACT_ARITHMETIC_POLICY<T,EXACT_RATIONAL<T> >
{
    typedef EXACT_RATIONAL<T> SUM_TYPE;
    typedef EXACT_RATIONAL<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
template<class T>
struct EXACT_ARITHMETIC_POLICY<EXACT_FLOAT<T>,EXACT_RATIONAL<T> >
{
    typedef EXACT_RATIONAL<T> SUM_TYPE;
    typedef EXACT_RATIONAL<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
template<class T>
struct EXACT_ARITHMETIC_POLICY<EXACT_RATIONAL<T>,EXACT_RATIONAL<T> >
{
    typedef EXACT_RATIONAL<T> SUM_TYPE;
    typedef EXACT_RATIONAL<T> PRODUCT_TYPE;
    typedef EXACT_RATIONAL<T> QUOTIENT_TYPE;
};
//#####################################################################
}
#endif
