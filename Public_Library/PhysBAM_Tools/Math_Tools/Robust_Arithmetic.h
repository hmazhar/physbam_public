//#####################################################################
// Copyright 2006-2007, Don Hatch, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __safe_arithmetic__
#define __safe_arithmetic__

#include <PhysBAM_Tools/Math_Tools/exchange_sort.h>
#include <cfloat>
#include <cmath>
namespace PhysBAM{

using ::std::abs;

template<class T>
inline T Robust_Multiply(const T a,const T b)
{T abs_a=abs(a),abs_b=abs(b);
if(abs_b<1 || abs_a<FLT_MAX/abs_b) return a*b;
else return (a>0)==(b>=0)?FLT_MAX:-FLT_MAX;}

template<class T>
inline T Robust_Divide(const T a,const T b)
{T abs_a=abs(a),abs_b=abs(b);
if(abs_b==FLT_MAX) return T();
if(abs_b>1 || abs_a<FLT_MAX*abs_b) return a/b;
else return (a>0)==(b>=0)?FLT_MAX:-FLT_MAX;}

template<class T>
inline T Robust_Inverse(const T a)
{T abs_a=abs(a);
if(abs_a==FLT_MAX) return T();
if(abs_a>1 || FLT_MAX*abs_a>1) return 1/a;
else return (a>=0)?FLT_MAX:-FLT_MAX;}

template<class T>
inline T Pseudo_Inverse(const T a)
{if(a==0) return T();else return (T)1/a;}

template<class T1,class T2>
inline T1 Pseudo_Divide(const T1& a,const T2& b)
{if(b==0) return T1();else return a/b;}

template<class T>
inline T Robust_Harmonic_Mean(T a,T b)
{exchange_sort(a,b);assert(a>=0);
return b>0?a/(1+a/b)*2:0;}

}
#endif
