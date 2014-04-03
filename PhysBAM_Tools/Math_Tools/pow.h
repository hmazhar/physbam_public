//#####################################################################
// Copyright 2006, Geoffrey Irving, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function pow
//#####################################################################
//
// pow is slow
//
//#####################################################################
#ifndef __pow__
#define __pow__

#include <PhysBAM_Tools/Math_Tools/cbrt.h>
namespace PhysBAM{

using ::std::sqrt;
using ::std::pow;

template<class T,int numerator,int denominator=1> struct POW_HELPER;

template<class T> struct POW_HELPER<T,-3>{static T pow(const T a){return 1/(a*a*a);}};
template<class T> struct POW_HELPER<T,-2>{static T pow(const T a){return 1/(a*a);}};
template<class T> struct POW_HELPER<T,-1>{static T pow(const T a){return 1/a;}};
template<class T> struct POW_HELPER<T,0>{static T pow(const T a){return 1;}};
template<class T> struct POW_HELPER<T,1,2>{static T pow(const T a){return sqrt(a);}};
template<class T> struct POW_HELPER<T,1,3>{static T pow(const T a){return cbrt(a);}};
template<class T> struct POW_HELPER<T,1>{static T pow(const T a){return a;}};
template<class T> struct POW_HELPER<T,2>{static T pow(const T a){return a*a;}};
template<class T> struct POW_HELPER<T,3>{static T pow(const T a){return a*a*a;}};

template<int numerator,class T> inline T pow(const T a)
{
    return POW_HELPER<T,numerator,1>::pow(a);
}

template<int numerator,int denominator,class T> inline T pow(const T a)
{
    return POW_HELPER<T,numerator,denominator>::pow(a);
}

}
#endif
