//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ROBUST_FUNCTIONS
//#####################################################################
#ifndef __ROBUST_FUNCTIONS__
#define __ROBUST_FUNCTIONS__

#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <cmath>
namespace PhysBAM{

using ::std::abs;
using ::std::atan2;
using ::std::sin;

template<class T>
inline T sinc(const T x) // sin(x)/x
{if(abs(x)<1e-8) return 1;return sin(x)/x;}

template<class T>
inline T one_minus_cos_x_over_x_squared(const T x) // (1-cos(x))/x^2
{return (T).5*sqr(sinc((T).5*x));}

template<class T>
inline T one_minus_cos_x_over_x(const T x) // (1-cos(x))/x
{return one_minus_cos_x_over_x_squared(x)*x;}

template<class T>
inline T atan2_y_x_over_y(const T y,const T x) // atan2(y,x)/y
{if(abs(y)<1e-8) return 1;return atan2(y,x)/y;}
}
#endif
