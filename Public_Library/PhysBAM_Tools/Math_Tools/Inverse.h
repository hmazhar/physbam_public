//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Inverse
//#####################################################################
#ifndef __Inverse__
#define __Inverse__

#include <cassert>
#include <cfloat>
#include <cmath>
namespace PhysBAM{

using ::std::abs;

template<class T_MATRIX> inline T_MATRIX
Inverse(const T_MATRIX& x)
{
    return x.Inverse();
}

inline float
Inverse(const float x)
{assert(abs(x)>=FLT_MIN);return 1/x;}

inline double
Inverse(const double x)
{assert(abs(x)>=DBL_MIN);return 1/x;}

struct INT_INVERSE
{
    int a;
};

inline INT_INVERSE
Inverse(const int x)
{assert(x!=0);INT_INVERSE r;r.a=x;return r;}

inline int
operator*(const int x,const INT_INVERSE y)
{return x/y.a;}

inline int&
operator*=(int& x,const INT_INVERSE y)
{return x/=y.a;}

}
#endif

