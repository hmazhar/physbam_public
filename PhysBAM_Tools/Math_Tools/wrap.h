//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Igor Neverov, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function wrap
//#####################################################################
//
// wrap(i,n) adds a multiple of n to i to bring it into the set [0;n)
// i + k*n is in [0,n)
//               
//#####################################################################
#ifndef __wrap__
#define __wrap__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <cmath>
namespace PhysBAM{

inline int wrap(const int i,const int n)
{if(i>=0) return i%n;
else return n-(-i)%n;} // i < 0 case

inline float wrap(const float value,const float lower,const float upper)
{float r=fmod(value-lower,upper-lower);if(r<0) return r+upper;else return r+lower;}

inline double wrap(const double value,const double lower,const double upper)
{double r=fmod(value-lower,upper-lower);if(r<0) return r+upper;else return r+lower;}

template<class T,int d>
inline VECTOR<T,d> wrap(const VECTOR<T,d>& value,const VECTOR<T,d>& lower,const VECTOR<T,d>& upper)
{VECTOR<T,d> result;for(int i=1;i<=result.Size();i++) result(i)=wrap(value(i),lower(i),upper(i));return result;}

}
#endif
